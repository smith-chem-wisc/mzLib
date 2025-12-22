using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;
using Quantification.Strategies;
using System;
using MassSpectrometry;
using MzLibUtil;

namespace Quantification;

/// <summary>
/// A quantification engine that performs the following steps:
/// 0) Writes raw information snapshot (if enabled)
/// 1) Pivot: Covert the raw info (long, one row for every PSM) to a wide QuantMatrix 
/// 2) Normalize the Psm Matrix
/// 3) Roll up to peptides
///     3a) Map the PSMs to peptides, creating a Dictionary<IBioPolymerWithSetMods, List<int>> Map mapping peptides to the indices of their PSMs in the QuantMatrix
///     3b) Roll-up. The roll-up strategy will take in a QuantMatrix of PSMs and the > map, and output a Peptide QuantMatrix
/// 4) Normalize the peptide matrix
/// 5) 
/// 5) Roll up to proteins
///     5a) Map the peptides to proteins, creating a Dictionary<IBioPolymerGroup, List<int>> Map mapping proteins to the indices of their peptides in the QuantMatrix
///     5b) Roll-up. The roll-up strategy will take in a QuantMatrix of peptides and the map, and output a Protein QuantMatrix
/// 6) Normalize the protein matrix
/// normalize -> rollups -> optional writes.
/// Does not depend on MetaMorpheus; MetaMorpheus supplies the objects.
/// </summary>
public sealed class QuantificationEngine
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

        // 1) Write immutable raw snapshot
        if (Parameters.WriteRawInformation)
            QuantificationWriter.WriteRawData(SpectralMatches, Parameters.OutputDirectory);

        // 2) Pivot: Convert long format (one row per PSM) to wide format (QuantMatrix)
        var psmMatrix = Pivot(SpectralMatches, ExperimentalDesign);

        // 3) Normalize PSM matrix
        var psmMatrixNorm = Parameters.SpectralMatchNormalizationStrategy.NormalizeIntensities(psmMatrix);

        // 4) Roll up to peptides
        var peptideMap = GetPsmToPeptideMap(psmMatrixNorm);
        var peptideMatrix = Parameters.SpectralMatchToPeptideRollUpStrategy
            .RollUp(psmMatrixNorm, peptideMap);

        // 5) Normalize peptide matrix 
        var peptideMatrixNorm = Parameters.PeptideNormalizationStrategy.NormalizeIntensities(peptideMatrix);

        // 6) Roll up to protein groups
        var proteinMap = Parameters.UseSharedPeptidesForProteinQuant ?
            GetAllPeptideToProteinMap(peptideMatrixNorm) :
            GetUniquePeptideToProteinMap(peptideMatrixNorm);
        var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
            .RollUp(peptideMatrixNorm, proteinMap);

        // 7) Normalize protein group matrix
        var proteinMatrixNorm = Parameters.ProteinNormalizationStrategy.NormalizeIntensities(proteinMatrix);

        // 8) Write results (If enabled)
        if(Parameters.WritePeptideInformation)
            QuantificationWriter.WritePeptideMatrix(peptideMatrixNorm, Parameters.OutputDirectory);
        if(Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrixNorm, Parameters.OutputDirectory);

        return new QuantificationResults
        {
            Summary = "Quantification completed successfully."
        };
    }

    public SpectralMatchMatrix Pivot(List<ISpectralMatch> spectralMatches, IExperimentalDesign experimentalDesign)
    {
        var sampleInfos = GetOrderedSampleInfos(experimentalDesign, spectralMatches, out var filePathToArrayPositionDict);

        var quantifiedMatches = spectralMatches
            .Where(sm => sm.QuantValues != null)
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
                 intensities[arrayPosition] += spectralMatch.QuantValues[i];
            }
            // Copy the summed intensities to the matrix
            smMatrix.SetRow(spectralMatch, intensities);
        }
        return smMatrix;
    }

    /// <summary>
    /// Orders first by file name (alphabetically), then by sample info for each file
    /// </summary>
    /// <param name="experimentalDesign"></param>
    /// <param name="spectralMatches"></param>
    /// <param name="filePathToArrayPositionDict"> Dictionary mapping file paths to their corresponding zero-indexed array positions in the ordered sample list </param>
    /// <returns>A list of ISampleInfo objects. This list can be used designate columns in a QuantMatrix object </returns>
    public List<ISampleInfo> GetOrderedSampleInfos(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches,
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

    public Dictionary<IBioPolymerWithSetMods, List<int>> GetPsmToPeptideMap(QuantMatrix<ISpectralMatch> psmMatrix)
    {
        var peptideToPsmMap = new Dictionary<IBioPolymerWithSetMods, List<int>>();
        for (int i = 0; i < psmMatrix.RowKeys.Count; i++)
        {
            var sm = psmMatrix.RowKeys[i];
            var peptide = sm.GetIdentifiedBioPolymersWithSetMods().First(); // Assumes unambiguous mapping
            if (!peptideToPsmMap.ContainsKey(peptide))
            {
                peptideToPsmMap[peptide] = new List<int>();
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
    public Dictionary<IBioPolymerGroup, List<int>> GetUniquePeptideToProteinMap(QuantMatrix<IBioPolymerWithSetMods> peptideMatrix)
    {
        var proteinToPeptideMap = new Dictionary<IBioPolymerGroup, List<int>>();
        
        // Initialize empty lists for each protein group
        foreach (var protein in BioPolymerGroups)
        {
            proteinToPeptideMap[protein] = new List<int>();
        }

        // Create a dictionary that maps each unique peptide to its corresponding protein group
        // Each peptide belongs to exactly one protein group (though it may be shared across proteins within that group)
        var peptideToProteinMap = new Dictionary<IBioPolymerWithSetMods, IBioPolymerGroup>();
        foreach (var proteinGroup in BioPolymerGroups)
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

    // TODO: Implement this method to include all peptides (not just unique ones) in the mapping
    public Dictionary<IBioPolymerGroup, List<int>> GetAllPeptideToProteinMap(
        QuantMatrix<IBioPolymerWithSetMods> peptideMatrix)
    {
        throw new NotImplementedException();
    }
}
