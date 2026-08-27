using Omics;
using Omics.SpectralMatch;
using Omics.BioPolymerGroup;
using MassSpectrometry;
using MzLibUtil;

namespace Quantification;

/// <summary>
/// A quantification engine that performs the following main (numbered) and ancillary steps:
/// 1) Creates one SpectralMatchMatrix per file
///     1a) Normalize the SpectralMatch Matrix for each file
/// 2) Roll up to peptides for each file
///     2a) Map the PSMs to peptides, creating a Dictionary of IBioPolymerWithSetMods to List of int, mapping peptides to the indices of their PSMs in the QuantMatrix
///     2b) Roll-up. The roll-up strategy will take in a QuantMatrix of PSMs and the map, and output a Peptide QuantMatrix
///     2c) Combine the per-file peptide matrices into a single matrix spanning all files, with missing values filled in as 0s
/// 3) Normalize the peptide matrix
/// 4) Collapse the peptide matrix, combining fractions and technical replicates
/// * Writes peptide information (if enabled)
/// 5) Roll up to proteins
///     5a) Map the peptides to proteins, creating a Dictionary of IBioPolymerGroup to List of int, mapping proteins to the indices of their peptides in the QuantMatrix
///     5b) Roll-up. The roll-up strategy will take in a QuantMatrix of peptides and the map, and output a Protein QuantMatrix
/// 6) Normalize the protein matrix
/// * Writes protein information (if enabled)
/// </summary>
public class QuantificationEngine
{
    public QuantificationParameters Parameters { get; init; }
    public IExperimentalDesign ExperimentalDesign { get; init; }
    internal List<ISpectralMatch> SpectralMatches { get; init; }
    internal List<IBioPolymerWithSetMods> ModifiedBioPolymers { get; init; }
    internal List<IBioPolymerGroup> BioPolymerGroups { get; init; }
    
    public QuantificationEngine(
        QuantificationParameters parameters,
        IExperimentalDesign experimentalDesign,
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> modifiedBioPolymers,
        List<IBioPolymerGroup> bioPolymerGroups)
    {
        Parameters = parameters;
        ExperimentalDesign = experimentalDesign;
        SpectralMatches = spectralMatches;
        ModifiedBioPolymers = modifiedBioPolymers;
        BioPolymerGroups = bioPolymerGroups;
    }

    /// <summary>
    /// The directory this run writes to: <see cref="QuantificationParameters.OutputDirectory"/> when the
    /// caller set one, otherwise the directory the source data files live in. Null until
    /// <see cref="ValidateEngine"/> has run, and null for a run with all three write flags off.
    /// </summary>
    public string ResolvedOutputDirectory { get; private set; }

    public QuantificationResults Run()
    {
        return Run(out var proteinMatrix);
    }

    public QuantificationResults Run(out QuantMatrix<IBioPolymerGroup> proteinMatrix)
    {
        proteinMatrix = null;

        // 0) Validate engine state
        if (!ValidateEngine(out QuantificationResults badResults))
        {
            return badResults;
        }

        var writtenFiles = new List<string>();

        // Write immutable raw snapshot
        Task<string> rawWriter = null;
        if (Parameters.WriteRawInformation)
        {
            rawWriter = Task.Run(() =>
                QuantificationWriter.WriteRawData(SpectralMatches, ResolvedOutputDirectory));
        }

        RunPeptideQuant(out var peptideMatrix);

        Task<string> peptideWriter = null;
        if (Parameters.WritePeptideInformation)
        {
            peptideWriter = Task.Run(() =>
                QuantificationWriter.WritePeptideMatrix(peptideMatrix, ResolvedOutputDirectory));
        }

        RunProteinQuant(peptideMatrix, out proteinMatrix);

        // 9) Write protein results (If enabled) No need to spin up a task, as this is the last step and we need to wait for it to complete before returning results anyway
        string proteinFile = null;
        if (Parameters.WriteProteinInformation)
            proteinFile = QuantificationWriter.WriteProteinGroupMatrix(proteinMatrix, ResolvedOutputDirectory);

        if (rawWriter != null) writtenFiles.Add(rawWriter.GetAwaiter().GetResult());
        if (peptideWriter != null) writtenFiles.Add(peptideWriter.GetAwaiter().GetResult());
        if (proteinFile != null) writtenFiles.Add(proteinFile);

        // 10) Deliver. Write the final values back onto every entity that opts in via
        //     IHasSampleIntensities, so that the existing BioPolymerGroup writers can render them,
        //     and collect the same values into the returned results object.
        OverwriteSampleIntensities(peptideMatrix);
        OverwriteSampleIntensities(proteinMatrix);

        // Report the matches the peptide map had to drop, so an unexpectedly small result set has a
        // stated cause rather than being silently smaller than the search.
        // Same predicate the peptide map applies, so the count cannot drift from what was dropped.
        int ambiguousExcluded = SpectralMatches
            .Count(sm => sm.Intensities != null && SingleIdentifiedBioPolymerOrNull(sm) == null);

        return new QuantificationResults
        {
            Summary = ambiguousExcluded == 0
                ? "Quantification completed successfully."
                : $"Quantification completed successfully. {ambiguousExcluded} spectral match(es) were " +
                  "excluded because they did not identify exactly one biopolymer.",
            Success = true,
            AmbiguousSpectralMatchesExcluded = ambiguousExcluded,
            Samples = proteinMatrix.ColumnKeys.ToList(),
            PeptideIntensities = BuildIntensityTable(peptideMatrix),
            ProteinIntensities = BuildIntensityTable(proteinMatrix),
            WrittenFiles = writtenFiles,
            OutputDirectory = ResolvedOutputDirectory
        };
    }

    /// <summary>
    /// True when any of the three write flags is set, and the run therefore needs an output directory.
    /// </summary>
    private bool WritesAnything() =>
        Parameters.WriteRawInformation || Parameters.WritePeptideInformation || Parameters.WriteProteinInformation;

    /// <summary>
    /// Replaces the values on any row key that implements <see cref="IHasSampleIntensities"/> with the
    /// final matrix values, setting both <see cref="IHasSampleIntensities.SamplesForQuantification"/>
    /// (the column order) and <see cref="IHasSampleIntensities.IntensitiesBySample"/> (the values).
    /// </summary>
    /// <remarks>
    /// This is destructive by design — hence "Overwrite". Anything already on the entity is discarded,
    /// so an entity carrying values from an earlier run comes out holding this run's values only.
    /// The pre-normalization PSM values are not lost: <see cref="QuantificationWriter.WriteRawData"/>
    /// records them before any normalization or roll-up runs, which is what makes re-processing with
    /// different strategies possible without re-searching.
    ///
    /// Row keys that do not implement the interface are skipped, so this is safe to call on a peptide
    /// matrix even though <see cref="IBioPolymerWithSetMods"/> does not carry quantification state.
    /// Zero-valued cells are omitted: the matrix uses 0 to mean "not observed in this sample", and
    /// omitting them keeps <c>SampleGroupResult.HasIntensityData</c> meaningful.
    /// </remarks>
    internal static void OverwriteSampleIntensities<T>(QuantMatrix<T> matrix) where T : IEquatable<T>
    {
        if (matrix == null) return;

        for (int row = 0; row < matrix.RowCount; row++)
        {
            if (matrix.RowKeys[row] is not IHasSampleIntensities entity)
                continue;

            var intensities = new Dictionary<ISampleInfo, double>(matrix.ColumnCount);
            for (int col = 0; col < matrix.ColumnCount; col++)
            {
                double value = matrix.Matrix[row, col];
                if (value != 0)
                    intensities[matrix.ColumnKeys[col]] = value;
            }

            // Assign intensities before samples: BioPolymerGroup's setters invalidate the cached
            // SampleGroupResults, and we want the invalidation to be the last thing that happens.
            entity.IntensitiesBySample = intensities;
            entity.SamplesForQuantification = matrix.ColumnKeys.ToList();
        }
    }

    /// <summary>
    /// Projects a matrix into a row-key → (sample → intensity) table for <see cref="QuantificationResults"/>.
    /// Zero-valued cells are omitted, matching <see cref="OverwriteSampleIntensities{T}"/>.
    /// </summary>
    internal static Dictionary<T, Dictionary<ISampleInfo, double>> BuildIntensityTable<T>(QuantMatrix<T> matrix)
        where T : IEquatable<T>
    {
        var table = new Dictionary<T, Dictionary<ISampleInfo, double>>();
        if (matrix == null) return table;

        for (int row = 0; row < matrix.RowCount; row++)
        {
            var intensities = new Dictionary<ISampleInfo, double>(matrix.ColumnCount);
            for (int col = 0; col < matrix.ColumnCount; col++)
            {
                double value = matrix.Matrix[row, col];
                if (value != 0)
                    intensities[matrix.ColumnKeys[col]] = value;
            }
            table[matrix.RowKeys[row]] = intensities;
        }

        return table;
    }

    /// <summary>
    /// Checks Engine state for validity before running quantification.
    /// </summary>
    /// <param name="badResults"> Return quant results with descriptive Summary if problem was encountered; null otherwise</param>
    /// <returns> True if engine can be ran succesfully, false otherwise </returns>
    internal bool ValidateEngine(out QuantificationResults badResults)
    {
        badResults = null;
        if (ExperimentalDesign == null)
        {
            badResults = QuantificationResults.Failure("QuantificationEngine Error: Experimental design is null.");
            return false;
        }
        if(SpectralMatches.IsNullOrEmpty())
        {
            badResults = QuantificationResults.Failure("QuantificationEngine Error: No spectral matches provided for quantification.");
            return false;
        }
        if(ModifiedBioPolymers.IsNullOrEmpty())
        {
            badResults = QuantificationResults.Failure("QuantificationEngine Error: No modified biopolymers (peptides) provided for quantification.");
            return false;
        }
        if(BioPolymerGroups.IsNullOrEmpty())
        {
            badResults = QuantificationResults.Failure("QuantificationEngine Error: No biopolymer groups (proteins) provided for quantification.");
            return false;
        }
        // Writing is on by default, so a caller who never set OutputDirectory would otherwise scatter
        // files into the working directory. Fall back to the directory the data came from, and only
        // if even that cannot be worked out, say so rather than guessing.
        ResolvedOutputDirectory = null;
        if (WritesAnything())
        {
            if (!string.IsNullOrWhiteSpace(Parameters.OutputDirectory))
            {
                ResolvedOutputDirectory = Parameters.OutputDirectory;
            }
            else if (TryGetSourceFileDirectory(SpectralMatches, out string beside))
            {
                ResolvedOutputDirectory = beside;
            }
            else
            {
                badResults = QuantificationResults.Failure(
                    "QuantificationEngine Error: OutputDirectory is not set and no default could be " +
                    "derived, because the spectral matches' file paths are missing, relative, or span " +
                    "unrelated directories. Set OutputDirectory, or turn WriteRawInformation, " +
                    "WritePeptideInformation and WriteProteinInformation off to quantify without writing files.");
                return false;
            }
        }
        return true;
    }

    /// <summary>
    /// The default output directory: where the source data files live, so that output lands beside the
    /// data it came from when the caller did not name a directory.
    /// </summary>
    /// <remarks>
    /// Only absolute source paths are considered. A bare file name or a relative path resolves against
    /// the process's working directory, which is precisely the ambiguity this default exists to avoid,
    /// so such matches are ignored rather than silently anchored to wherever the caller happened to be.
    ///
    /// When the files span several directories the nearest common ancestor is used, so one search over
    /// <c>data/fraction1</c> and <c>data/fraction2</c> writes to <c>data</c>. A common ancestor that is
    /// a drive or filesystem root is rejected: files on unrelated volumes, or in unrelated top-level
    /// folders, have no shared home and the caller has to say where output goes.
    /// </remarks>
    /// <returns>True if a directory could be derived, in which case <paramref name="directory"/> holds it.</returns>
    internal static bool TryGetSourceFileDirectory(IEnumerable<ISpectralMatch> spectralMatches, out string directory)
    {
        directory = null;
        if (spectralMatches == null) return false;

        var directories = new List<string>();
        foreach (var match in spectralMatches)
        {
            string path = match?.FullFilePath;
            if (string.IsNullOrWhiteSpace(path)) continue;

            string dir;
            try
            {
                dir = Path.GetDirectoryName(path);
            }
            catch (ArgumentException)
            {
                continue; // malformed path; it contributes nothing to the default
            }

            if (string.IsNullOrEmpty(dir) || !Path.IsPathRooted(dir)) continue;
            if (!directories.Contains(dir, StringComparer.OrdinalIgnoreCase))
                directories.Add(dir);
        }

        if (directories.Count == 0) return false;
        if (directories.Count == 1)
        {
            directory = directories[0];
            return true;
        }

        string common = directories[0];
        for (int i = 1; i < directories.Count && common != null; i++)
            common = CommonAncestor(common, directories[i]);

        // Path.GetDirectoryName returns null only for a root, which is not a meaningful "beside the data".
        if (common == null || Path.GetDirectoryName(common) == null) return false;

        directory = common;
        return true;
    }

    /// <summary>
    /// The deepest directory that is a prefix of both paths, or null if they share nothing above the
    /// root. Comparison is per path segment, so <c>C:\dataset</c> and <c>C:\dataset2</c> share
    /// <c>C:\</c> rather than the character prefix <c>C:\dataset</c>.
    /// </summary>
    private static string CommonAncestor(string first, string second)
    {
        var separators = new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar };
        var a = first.Split(separators, StringSplitOptions.None);
        var b = second.Split(separators, StringSplitOptions.None);

        int shared = 0;
        while (shared < a.Length && shared < b.Length &&
               string.Equals(a[shared], b[shared], StringComparison.OrdinalIgnoreCase))
        {
            shared++;
        }

        if (shared == 0) return null; // different drives, or one rooted and one not

        // One shared segment is the root itself: "C:" from "C:\x", or "" from "/x".
        string joined = string.Join(Path.DirectorySeparatorChar.ToString(), a, 0, shared);
        return shared == 1 ? joined + Path.DirectorySeparatorChar : joined;
    }

    internal void RunPeptideQuant(out QuantMatrix<IBioPolymerWithSetMods> peptideMatrixNorm)
    {
        // 1) PivotByFile - one matrix per file
        var perFileMatrices = PivotByFile(SpectralMatches, ExperimentalDesign);

        // 2) Per-file PSM normalization
        var perFileNormalized = new Dictionary<string, QuantMatrix<ISpectralMatch>>();
        foreach (var kvp in perFileMatrices)
        {
            perFileNormalized[kvp.Key] = Parameters.SpectralMatchNormalizationStrategy
                .NormalizeIntensities(kvp.Value);
        }

        // 3) Per-file roll-up to peptides
        var perFilePeptideMatrices = new Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>>();
        foreach (var kvp in perFileNormalized)
        {
            var peptideMap = GetPsmToPeptideMap(kvp.Value, ModifiedBioPolymers);
            perFilePeptideMatrices[kvp.Key] = Parameters.SpectralMatchToPeptideRollUpStrategy
                .RollUp(kvp.Value, peptideMap);
        }

        // 4) Combine per-file peptide matrices
        var combinedPeptideMatrix = CombinePeptideMatrices(perFilePeptideMatrices, ExperimentalDesign);

        // 5) Normalize combined peptide matrix
        peptideMatrixNorm = Parameters.PeptideNormalizationStrategy
            .NormalizeIntensities(combinedPeptideMatrix);
    }

    internal void RunProteinQuant(QuantMatrix<IBioPolymerWithSetMods> peptideMatrixNorm,
        out QuantMatrix<IBioPolymerGroup> proteinMatrixNorm)
    {
        // 6) Collapse samples (technical replicates, fractions)
        peptideMatrixNorm = Parameters.CollapseStrategy.CollapseSamples(peptideMatrixNorm, Parameters.CollapseAggregationStrategy);

        // 7) Roll up to proteins
        var proteinMap = Parameters.UseSharedPeptidesForProteinQuant
            ? GetAllPeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups)
            : GetUniquePeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups);

        var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
            .RollUp(peptideMatrixNorm, proteinMap);

        // 8) Normalize protein matrix
        proteinMatrixNorm = Parameters.ProteinNormalizationStrategy
            .NormalizeIntensities(proteinMatrix);
    }

    /// <summary>
    /// Creates one SpectralMatchMatrix per file.
    /// Each matrix has rows = PSMs from that file, columns = channels within that file.
    /// This produces dense matrices (no sparse zeros) suitable for within-file normalization.
    /// For LFQ, each matrix will only have one column
    /// </summary>
    /// <param name="spectralMatches">All spectral matches across all files</param>
    /// <param name="experimentalDesign">Maps file names to channel ISampleInfo arrays</param>
    /// <returns>Dictionary mapping file path to its SpectralMatchMatrix</returns>
    public static Dictionary<string, SpectralMatchMatrix> PivotByFile(
        List<ISpectralMatch> spectralMatches,
        IExperimentalDesign experimentalDesign)
    {
        var result = new Dictionary<string, SpectralMatchMatrix>();

        // Filter to spectral matches with non-null Intensities and group by file path
        var quantified = spectralMatches
            .Where(sm => sm.Intensities != null)
            .GroupBy(sm => sm.FullFilePath)
            .OrderBy(g => g.Key);

        foreach (var fileGroup in quantified)
        {
            string filePath = fileGroup.Key;
            string fileName = Path.GetFileName(filePath);

            // Lookup the channel sample infos for this file
            if (!experimentalDesign.FileNameSampleInfoDictionary.TryGetValue(fileName, out var sampleInfoArray))
            {
                throw new KeyNotFoundException(
                    $"File name '{fileName}' not found in experimental design.");
            }

            // PSMs ordered by FullSequence for determinism
            var filePsms = fileGroup.OrderBy(sm => sm.FullSequence).ToList();

            var smMatrix = new SpectralMatchMatrix(filePsms, sampleInfoArray, experimentalDesign);

            // Copy Intensities directly — positional mapping: Intensities[i] → column[i]
            foreach (var sm in filePsms)
            {
                smMatrix.SetRow(sm, sm.Intensities);
            }

            result[filePath] = smMatrix;
        }

        return result;
    }

    /// <summary>
    /// Combines per-file peptide matrices into a single matrix spanning all files.
    /// Rows = union of all peptides across files.
    /// Columns = all channels from all files, ordered by file path then channel.
    /// Values = peptide intensity in that channel, or 0 if the peptide was not observed in that file.
    /// </summary>
    /// <param name="perFilePeptideMatrices">Dictionary of file path → peptide matrix for that file</param>
    /// <param name="experimentalDesign">The experimental design</param>
    /// <returns>A single PeptideMatrix covering all files and channels</returns>
    public static PeptideMatrix CombinePeptideMatrices(
        Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>> perFilePeptideMatrices,
        IExperimentalDesign experimentalDesign)
    {
        // 1. Collect all unique peptides across all per-file matrices, ordered for determinism
        var allPeptides = perFilePeptideMatrices.Values
            .SelectMany(m => m.RowKeys)
            .Distinct()
            .OrderBy(p => p.FullSequence)
            .ToList();

        // 2. Collect all column keys ordered by file path (alphabetically), then channel within each file
        var sortedFilePaths = perFilePeptideMatrices.Keys.OrderBy(fp => fp).ToList();
        var allColumnKeys = new List<ISampleInfo>();
        var fileColumnOffsets = new Dictionary<string, int>();
        foreach (var filePath in sortedFilePaths)
        {
            fileColumnOffsets[filePath] = allColumnKeys.Count;
            allColumnKeys.AddRange(perFilePeptideMatrices[filePath].ColumnKeys);
        }

        // 3. Create the combined PeptideMatrix
        var combined = new PeptideMatrix(allPeptides, allColumnKeys, experimentalDesign);

        // Pre-build peptide index map for O(1) lookup instead of O(n) IndexOf
        var peptideIndexMap = new Dictionary<IBioPolymerWithSetMods, int>(allPeptides.Count);
        for (int i = 0; i < allPeptides.Count; i++)
            peptideIndexMap[allPeptides[i]] = i;

        // 4. For each per-file matrix, copy values directly into the combined matrix at the correct offset
        foreach (var filePath in sortedFilePaths)
        {
            var fileMatrix = perFilePeptideMatrices[filePath];
            int colOffset = fileColumnOffsets[filePath];
            int numFileCols = fileMatrix.ColumnKeys.Count;

            for (int fileRow = 0; fileRow < fileMatrix.RowKeys.Count; fileRow++)
            {
                var peptide = fileMatrix.RowKeys[fileRow];
                if (!peptideIndexMap.TryGetValue(peptide, out int combinedRow))
                    continue;

                // Copy directly from source matrix to combined matrix — no GetRow/SetRow allocation
                for (int col = 0; col < numFileCols; col++)
                {
                    combined.Matrix[combinedRow, colOffset + col] = fileMatrix.Matrix[fileRow, col];
                }
            }
        }

        return combined;
    }

    /// <summary>
    /// Creates a mapping from each specified modified biopolymer to a list of indices that identify the position of corresponding
    /// spectral matches in the smMatrix
    /// </summary>
    /// <remarks>A spectral match is quantified only when it identifies exactly one modified biopolymer.
    /// An ambiguous match -- one that identifies several -- is excluded rather than attributed to whichever
    /// biopolymer happens to be enumerated first, and a match that identifies none is excluded rather than
    /// throwing. Biopolymers not present in the input list are ignored.
    /// <see cref="QuantificationResults.AmbiguousSpectralMatchesExcluded"/> reports how many were dropped.</remarks>
    /// <param name="smMatrix">The matrix containing spectrum matches to be mapped to their corresponding modified bioPolymer.</param>
    /// <param name="modifiedBioPolymers">The list of modified bioPolymers for which to generate the mapping.
    /// Only SMs corresponding to these bioPolymers are included in the result.</param>
    /// <returns>A dictionary mapping each modified bioPolymer in the input list to a list of indices of PSMs in the matrix that
    /// are associated with it. If a bioPolymer has no corresponding PSMs, its list will be empty.</returns>
    /// <summary>
    /// The single biopolymer a spectral match identifies, or null when it identifies none or several.
    /// This is the unambiguous filter, in one place, so that the quantified set and the count reported
    /// as <see cref="QuantificationResults.AmbiguousSpectralMatchesExcluded"/> cannot disagree.
    /// </summary>
    /// <remarks>
    /// Distinct, not raw count. A match that names the same biopolymer more than once -- once per
    /// protein it maps to, for instance -- identifies one peptide in substance, and dropping it as
    /// ambiguous would discard a perfectly good measurement. Nulls are ignored for the same reason.
    /// Take(2) still short-circuits, so a long ambiguity list is not enumerated.
    /// </remarks>
    internal static IBioPolymerWithSetMods SingleIdentifiedBioPolymerOrNull(ISpectralMatch spectralMatch)
    {
        var identified = spectralMatch.GetIdentifiedBioPolymersWithSetMods();
        if (identified == null)
        {
            return null;
        }

        var distinct = identified
            .Where(bp => bp != null)
            .Distinct()
            .Take(2)
            .ToList();

        return distinct.Count == 1 ? distinct[0] : null;
    }

    public static Dictionary<IBioPolymerWithSetMods, List<int>> GetPsmToPeptideMap(QuantMatrix<ISpectralMatch> smMatrix, List<IBioPolymerWithSetMods> modifiedBioPolymers)
    {
        var peptideToPsmMap = new Dictionary<IBioPolymerWithSetMods, List<int>>();
        foreach (var bp in modifiedBioPolymers)
        {
            peptideToPsmMap[bp] = new List<int>();
        }
        for (int i = 0; i < smMatrix.RowKeys.Count; i++)
        {
            var sm = smMatrix.RowKeys[i];

            // Only unambiguous matches are quantified.
            var peptide = SingleIdentifiedBioPolymerOrNull(sm);
            if (peptide == null)
            {
                continue;
            }

            if (!peptideToPsmMap.ContainsKey(peptide))
            {
                continue;
            }
            peptideToPsmMap[peptide].Add(i);
        }
        return peptideToPsmMap;
    }

    /// <summary>
    /// Creates a mapping from each protein group to the list of row indices in the peptide matrix that correspond to
    /// peptides uniquely assigned to that group.
    /// </summary>
    /// <remarks>Each peptide is assigned to exactly one protein group, even if it is shared among multiple
    /// proteins within that group. The returned mapping can be used to efficiently retrieve all peptides associated
    /// with a given protein group from the peptide matrix.</remarks>
    /// <param name="peptideMatrix">A matrix containing peptides as row keys, where each peptide is associated with a protein group.</param>
    /// <returns>A dictionary that maps each protein group to a list of integer indices. Each list contains the row indices in
    /// the peptide matrix for peptides uniquely assigned to the corresponding protein group. If a protein group has no
    /// assigned peptides, its list will be empty.</returns>
    public static Dictionary<IBioPolymerGroup, List<int>> GetUniquePeptideToProteinMap(QuantMatrix<IBioPolymerWithSetMods> peptideMatrix, List<IBioPolymerGroup> bioPolymerGroups)
    {
        var proteinToPeptideMap = new Dictionary<IBioPolymerGroup, List<int>>();
        
        // Initialize empty lists for each protein group
        foreach (var protein in bioPolymerGroups)
        {
            proteinToPeptideMap[protein] = new List<int>();
        }

        // Create a dictionary that maps each unique peptide to its corresponding protein group
        // Each peptide belongs to exactly one protein group (though it may be shared across proteins within that group)
        var peptideToProteinMap = new Dictionary<IBioPolymerWithSetMods, IBioPolymerGroup>();
        foreach (var proteinGroup in bioPolymerGroups)
        {
            foreach (var peptide in proteinGroup.UniqueBioPolymersWithSetMods)
            {
                peptideToProteinMap[peptide] = proteinGroup;
            }
        }

        // Iterate through the peptide matrix row keys and add each row index to its corresponding protein's list
        for (int i = 0; i < peptideMatrix.RowKeys.Count; i++)
        {
            var peptide = peptideMatrix.RowKeys[i];
            
            // Find which protein group this peptide belongs to
            if (peptideToProteinMap.TryGetValue(peptide, out var proteinGroup))
            {
                proteinToPeptideMap[proteinGroup].Add(i);
            }
        }

        return proteinToPeptideMap;
    }

    /// <summary>
    /// Creates a mapping from each protein group to the list of row indices in the peptide matrix that correspond to
    /// every peptide assigned to that group -- shared as well as unique.
    /// </summary>
    /// <remarks>A shared peptide belongs to more than one protein group, so its row index appears in more than one
    /// list. That is the difference from <see cref="GetUniquePeptideToProteinMap"/>, where each index appears once:
    /// a shared peptide's intensity contributes to every group it was assigned to. Indices are sorted, because
    /// <see cref="IBioPolymerGroup.AllBioPolymersWithSetMods"/> is a HashSet and its enumeration order is not
    /// guaranteed stable -- an unsorted list would make roll-up results depend on set ordering.</remarks>
    /// <param name="peptideMatrix">A matrix containing peptides as row keys.</param>
    /// <param name="bioPolymerGroups">The protein groups to map. Groups with no peptide in the matrix get an empty list.</param>
    /// <returns>A dictionary that maps each protein group to the sorted row indices of all of its peptides.</returns>
    public static Dictionary<IBioPolymerGroup, List<int>> GetAllPeptideToProteinMap(
        QuantMatrix<IBioPolymerWithSetMods> peptideMatrix, List<IBioPolymerGroup> bioPolymerGroups)
    {
        var proteinToPeptideMap = new Dictionary<IBioPolymerGroup, List<int>>();

        // Initialize empty lists for each protein group
        foreach (var protein in bioPolymerGroups)
        {
            proteinToPeptideMap[protein] = new List<int>();
        }

        // Index the matrix rows once, rather than scanning it per protein group
        var rowIndexByPeptide = new Dictionary<IBioPolymerWithSetMods, int>();
        for (int i = 0; i < peptideMatrix.RowKeys.Count; i++)
        {
            rowIndexByPeptide[peptideMatrix.RowKeys[i]] = i;
        }

        foreach (var proteinGroup in bioPolymerGroups)
        {
            var rowIndices = proteinToPeptideMap[proteinGroup];

            foreach (var peptide in proteinGroup.AllBioPolymersWithSetMods)
            {
                // A peptide the caller did not pass to the engine has no row to contribute
                if (rowIndexByPeptide.TryGetValue(peptide, out int rowIndex))
                {
                    rowIndices.Add(rowIndex);
                }
            }

            rowIndices.Sort();
        }

        return proteinToPeptideMap;
    }
}
