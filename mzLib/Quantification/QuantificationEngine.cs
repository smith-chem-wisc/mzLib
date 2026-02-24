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
///     2a) Map the PSMs to peptides, creating a Dictionary<IBioPolymerWithSetMods, List<int>> Map mapping peptides to the indices of their PSMs in the QuantMatrix
///     2b) Roll-up. The roll-up strategy will take in a QuantMatrix of PSMs and the > map, and output a Peptide QuantMatrix
///     2c) Combine the per-file peptide matrices into a single matrix spanning all files, with missing values filled in as 0s
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

        // Write immutable raw snapshot 
        Task rawWriter = null;
        if (Parameters.WriteRawInformation)
        {
            rawWriter = Task.Run(async () =>
            {
                QuantificationWriter.WriteRawData(SpectralMatches, Parameters.OutputDirectory);
            });
        }

        RunPeptideQuant(out var peptideMatrix);

        Task peptideWriter = null;
        if (Parameters.WritePeptideInformation)
        {
            peptideWriter = Task.Run(async () =>
            {
                QuantificationWriter.WritePeptideMatrix(peptideMatrix, Parameters.OutputDirectory);
            });
        }

        RunProteinQuant(peptideMatrix, out proteinMatrix);

        // 9) Write protein results (If enabled) No need to spin up a task, as this is the last step and we need to wait for it to complete before returning results anyway
        if (Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrix, Parameters.OutputDirectory);

        if (rawWriter != null) Task.WaitAll(rawWriter);
        if (peptideWriter != null) Task.WaitAll(peptideWriter);
        return new QuantificationResults
        {
            Summary = "Quantification completed successfully."
        };
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
