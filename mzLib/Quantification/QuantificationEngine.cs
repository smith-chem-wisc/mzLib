using Omics;
using Omics.SpectralMatch;
using Omics.BioPolymerGroup;
using MassSpectrometry;
using MzLibUtil;

namespace Quantification;

/// <summary>
/// A quantification engine that performs the following main (numbered) and ancillary steps:
/// * Pivot: Covert the raw info (long, one row for every PSM) to a wide QuantMatrix 
/// * Writes raw information snapshot (if enabled)
/// 1) Normalize the SpectralMatch Matrix
/// 2) Roll up to peptides
///     2a) Map the PSMs to peptides, creating a Dictionary<IBioPolymerWithSetMods, List<int>> Map mapping peptides to the indices of their PSMs in the QuantMatrix
///     2b) Roll-up. The roll-up strategy will take in a QuantMatrix of PSMs and the > map, and output a Peptide QuantMatrix
/// 3) Normalize the peptide matrix
/// 4) Collapse the peptide matrix, combining fractions and technical replicates
/// * Writes peptide information (if enabled)
/// 5) Roll up to proteins
///     5a) Map the peptides to proteins, creating a Dictionary<IBioPolymerGroup, List<int>> Map mapping proteins to the indices of their peptides in the QuantMatrix
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
    /// Checks Engine state for validity before running quantification.
    /// </summary>
    /// <param name="badResults"> Return quant results with descriptive Summary if problem was encountered; null otherwise</param>
    /// <returns> True if engine can be ran succesfully, false otherwise </returns>
    internal bool ValidateEngine(out QuantificationResults badResults)
    {
        badResults = null;
        if (ExperimentalDesign == null)
        {
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: Experimental design is null."
            };
            return false;
        }
        if(SpectralMatches.IsNullOrEmpty())
        {
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: No spectral matches provided for quantification."
            };
            return false;
        }
        if(ModifiedBioPolymers.IsNullOrEmpty())
        {
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: No modified biopolymers (peptides) provided for quantification."
            };
            return false;
        }
        if(BioPolymerGroups.IsNullOrEmpty())
        {
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: No biopolymer groups (proteins) provided for quantification."
            };
            return false;
        }
        return true;
    }
    
    public QuantificationResults Run()
    {
        // 0) Validate engine state
        if(!ValidateEngine(out QuantificationResults badResults))
        {
            return badResults;
        }

        // Write immutable raw snapshot 
        Task rawWriter = null;
        if (Parameters.WriteRawInformation)
        {
            rawWriter = Task.Run(async () =>
            {
                QuantificationWriter.WriteRawData(SpectralMatches, Parameters.OutputDirectory);
            });
        }

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


        // 3) Normalize peptide matrix 
        var peptideMatrixNorm = Parameters.PeptideNormalizationStrategy.NormalizeIntensities(combinedPeptideMatrix);

        // 4) Collapse the peptide matrix to combine fractions and technical replicates
        peptideMatrixNorm = Parameters.CollapseStrategy.CollapseSamples(peptideMatrixNorm);

        // Write peptide results (If enabled)
        if(Parameters.WritePeptideInformation) // TODO: Make this happen asynchronously
            QuantificationWriter.WritePeptideMatrix(peptideMatrixNorm, Parameters.OutputDirectory);

        // 5) Roll up to protein groups
        var proteinMap = Parameters.UseSharedPeptidesForProteinQuant ?
            GetAllPeptideToProteinMap(peptideMatrixNorm) :
            GetUniquePeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups);
        var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
            .RollUp(peptideMatrixNorm, proteinMap);

        // 6) Normalize protein group matrix
        var proteinMatrixNorm = Parameters.ProteinNormalizationStrategy.NormalizeIntensities(proteinMatrix);

        // 8) Write protein results (If enabled)
        if(Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrixNorm, Parameters.OutputDirectory);

        if (rawWriter != null) Task.WaitAll(rawWriter);
        return new QuantificationResults
        {
            Summary = "Quantification completed successfully."
        };
    }

    /// <summary>
    /// Runs the TMT-specific quantification pipeline.
    /// Processes each file independently (pivot, normalize, roll-up) before combining.
    /// </summary>
    public QuantificationResults RunTmt()
    {
        if (!ValidateEngine(out QuantificationResults badResults))
            return badResults;

        RunTmtCore(out _);
        return new QuantificationResults
        {
            Summary = "TMT Quantification completed successfully."
        };
    }

    /// <summary>
    /// This is only for testing until the writers are up and running. Returns the protein matrix
    /// </summary>
    /// <returns></returns>
    /// <exception cref="MzLibException"></exception>
    internal QuantMatrix<IBioPolymerGroup> RunAndReturnProteinMatrix()
    {
        // 0) Validate engine state
        if (!ValidateEngine(out QuantificationResults badResults))
        {
            throw new MzLibException("I don't know what went wrong here");
        }

        // Create a matrix from the spectral matches by converting from long format (one row per PSM) to wide format (QuantMatrix)
        var psmMatrix = Pivot(SpectralMatches, ExperimentalDesign);

        // Write immutable raw snapshot 
        // TODO: Make this happen asynchronously
        if (Parameters.WriteRawInformation)
            QuantificationWriter.WriteRawData(SpectralMatches, Parameters.OutputDirectory);

        // 1) Normalize PSM matrix
        var psmMatrixNorm = Parameters.SpectralMatchNormalizationStrategy.NormalizeIntensities(psmMatrix);

        // 2) Roll up to peptides
        var peptideMap = GetPsmToPeptideMap(psmMatrixNorm, ModifiedBioPolymers);
        var peptideMatrix = Parameters.SpectralMatchToPeptideRollUpStrategy
            .RollUp(psmMatrixNorm, peptideMap);

        // 3) Normalize peptide matrix 
        var peptideMatrixNorm = Parameters.PeptideNormalizationStrategy.NormalizeIntensities(peptideMatrix);

        // 4) Collapse the peptide matrix to combine fractions and technical replicates
        peptideMatrixNorm = Parameters.CollapseStrategy.CollapseSamples(peptideMatrixNorm);

        // Write peptide results (If enabled)
        if (Parameters.WritePeptideInformation) // TODO: Make this happen asynchronously
            QuantificationWriter.WritePeptideMatrix(peptideMatrixNorm, Parameters.OutputDirectory);

        // 5) Roll up to protein groups
        var proteinMap = Parameters.UseSharedPeptidesForProteinQuant ?
            GetAllPeptideToProteinMap(peptideMatrixNorm) :
            GetUniquePeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups);
        var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
            .RollUp(peptideMatrixNorm, proteinMap);

        // 6) Normalize protein group matrix
        var proteinMatrixNorm = Parameters.ProteinNormalizationStrategy.NormalizeIntensities(proteinMatrix);

        // 8) Write protein results (If enabled)
        if (Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrixNorm, Parameters.OutputDirectory);

        return proteinMatrixNorm;
    }


    /// <summary>
    /// Runs the TMT pipeline and returns the final protein matrix for testing/inspection.
    /// </summary>
    internal QuantMatrix<IBioPolymerGroup> RunTmtAndReturnProteinMatrix()
    {
        if (!ValidateEngine(out QuantificationResults _))
            throw new MzLibException("QuantificationEngine validation failed.");

        RunTmtCore(out var proteinMatrixNorm);
        return proteinMatrixNorm;
    }

    /// <summary>
    /// Core TMT pipeline: per-file normalization, roll-up, combine, normalize, collapse, protein roll-up.
    /// </summary>
    private void RunTmtCore(out QuantMatrix<IBioPolymerGroup> proteinMatrixNorm)
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
        var peptideMatrixNorm = Parameters.PeptideNormalizationStrategy
            .NormalizeIntensities(combinedPeptideMatrix);

        // 6) Collapse samples (technical replicates, fractions)
        peptideMatrixNorm = Parameters.CollapseStrategy.CollapseSamples(peptideMatrixNorm);

        // 7) Roll up to proteins
        var proteinMap = Parameters.UseSharedPeptidesForProteinQuant
            ? GetAllPeptideToProteinMap(peptideMatrixNorm)
            : GetUniquePeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups);
        var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
            .RollUp(peptideMatrixNorm, proteinMap);

        // 8) Normalize protein matrix
        proteinMatrixNorm = Parameters.ProteinNormalizationStrategy
            .NormalizeIntensities(proteinMatrix);
    }

    /// <summary>
    /// Creates a pivoted matrix of quantified spectral matches, organizing intensity values according to the specified
    /// experimental design.
    /// </summary>
    /// <remarks>Only spectral matches with non-null quantification values are included in the resulting
    /// matrix. The order of samples and the mapping of files to matrix columns are determined by the experimental
    /// design. This method is typically used to prepare data for downstream quantitative analysis.</remarks>
    /// <param name="spectralMatches">A list of spectral matches to be included in the matrix. Each match must
    /// have quantification values for inclusion.</param>
    /// <param name="experimentalDesign">The experimental design that defines how samples and files are mapped to matrix columns.</param>
    /// <returns>A SpectralMatchMatrix containing the quantified intensity values for each spectral match,
    /// arranged according to the experimental design.</returns>
    /// <exception cref="KeyNotFoundException">Thrown if a spectral match references a file path that is not present in the experimental design mapping.</exception>
    public static SpectralMatchMatrix Pivot(List<ISpectralMatch> spectralMatches, IExperimentalDesign experimentalDesign)
    {
        var sampleInfos = GetOrderedSampleInfos(experimentalDesign, spectralMatches, out var filePathToArrayPositionDict);

        var quantifiedMatches = spectralMatches
            .Where(sm => sm.Intensities != null)
            .OrderBy(sm => sm.FullSequence)
            .ToList();
        SpectralMatchMatrix smMatrix = new SpectralMatchMatrix(quantifiedMatches, sampleInfos, experimentalDesign);

        // create empty array to store intensities as they are summed, before they're copied to the matrix
        double[] intensities = new double[sampleInfos.Count];

        // Create a matrix where ever spectral match has its own row. 
        foreach (var spectralMatch in quantifiedMatches)
        {
            Array.Clear(intensities, 0, intensities.Length); // I'm pretty sure that values are copied when writing to the matrix,
                                                             // so this is fine. If it doesn't work, we can create a new array each time instead.

            // Each file is associated with one column (for LFQ) or multiple columns (for TMT) in the matrix
            // The arrayPositions dict enables mapping between file paths and their corresponding array positions in the matrix
            if (!filePathToArrayPositionDict.TryGetValue(spectralMatch.FullFilePath, out var arrayPositions))
            {
                throw new KeyNotFoundException($"File path '{spectralMatch.FullFilePath}' not found in file path to array position dictionary. This could indicate a problem with the ExperimentalDesign file");
            }

            // Copy the quant values to the appropriate position in the summedIntensities array
            for (var i = 0; i < arrayPositions.Count; i++)
            {
                 var arrayPosition = arrayPositions[i];
                 intensities[arrayPosition] = spectralMatch.Intensities[i];
            }
            // Copy the summed intensities to the matrix
            smMatrix.SetRow(spectralMatch, intensities);
        }
        return smMatrix;
    }

    /// <summary>
    /// Creates one SpectralMatchMatrix per file for TMT/isobaric data.
    /// Each matrix has rows = PSMs from that file, columns = channels within that file.
    /// This produces dense matrices (no sparse zeros) suitable for within-file normalization.
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
    /// Orders first by file name (alphabetically), then by sample info for each file
    /// </summary>
    /// <param name="experimentalDesign"></param>
    /// <param name="spectralMatches"></param>
    /// <param name="filePathToArrayPositionDict"> Dictionary mapping file paths to their corresponding zero-indexed array positions in the ordered sample list </param>
    /// <returns>A list of ISampleInfo objects. This list can be used designate columns in a QuantMatrix object </returns>
    public static List<ISampleInfo> GetOrderedSampleInfos(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches,
        out Dictionary<string, List<int>> filePathToArrayPositionDict)
    {
        // Start by ordering the file paths alphabetically        
        var filePaths = spectralMatches
            .Select(sm => sm.FullFilePath)
            .Distinct()
            .OrderBy(fp => fp)
            .ToList();

        List<ISampleInfo> orderedSamples = new();
        filePathToArrayPositionDict = new Dictionary<string, List<int>>();
        int zeroIndexedArrayPosition = 0;
        foreach (var filePath in filePaths)
        {
            var fileName = Path.GetFileName(filePath);
            filePathToArrayPositionDict[filePath] = new List<int>();
            if (experimentalDesign.FileNameSampleInfoDictionary.TryGetValue(fileName, out var sampleInfos))
            {
                orderedSamples.AddRange(sampleInfos);
                filePathToArrayPositionDict[filePath].AddRange(Enumerable.Range(zeroIndexedArrayPosition, sampleInfos.Length));
                zeroIndexedArrayPosition += sampleInfos.Length;
            }
            else
            {
                throw new KeyNotFoundException($"File name '{fileName}' not found in experimental design.");
            }
        }

        return orderedSamples;
    }

    /// <summary>
    /// Creates a mapping from each specified modified biopolymer to a list of indices that identify the position of corresponding
    /// spectral matches in the smMatrix
    /// </summary>
    /// <remarks>The mapping assumes that each PM in the matrix is associated with a single, unambiguous
    /// modified biopolymer. Biopolymers not present in the input list are ignored.</remarks>
    /// <param name="smMatrix">The matrix containing spectrum matches to be mapped to their corresponding modified bioPolymer.</param>
    /// <param name="modifiedBioPolymers">The list of modified bioPolymers for which to generate the mapping.
    /// Only SMs corresponding to these bioPolymers are included in the result.</param>
    /// <returns>A dictionary mapping each modified bioPolymer in the input list to a list of indices of PSMs in the matrix that
    /// are associated with it. If a bioPolymer has no corresponding PSMs, its list will be empty.</returns>
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
            var peptide = sm.GetIdentifiedBioPolymersWithSetMods().First(); // Assumes unambiguous mapping
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

    // TODO: Implement this method to that include all peptides (shared and unique) in the mapping
    public Dictionary<IBioPolymerGroup, List<int>> GetAllPeptideToProteinMap(
        QuantMatrix<IBioPolymerWithSetMods> peptideMatrix)
    {
        throw new NotImplementedException();
    }
}
