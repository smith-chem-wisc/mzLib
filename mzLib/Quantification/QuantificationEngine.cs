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
    public IRollUpStrategy RollUpStrategy { get; internal set; }
    public INormalizationStrategy NormalizationStrategy { get; internal set; }

    public QuantificationEngine(
        QuantificationParameters parameters,
        IExperimentalDesign experimentalDesign,
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> modifiedBioPolymers,
        List<IBioPolymerGroup> bioPolymerGroups)
    {
        Parameters = parameters;
        RollUpStrategy = Parameters.RollUpStrategy;
        NormalizationStrategy = Parameters.NormalizationStrategy;
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
        if(RollUpStrategy == null)
        {
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: Roll-up strategy is null."
            };
            return false;
        }
        if(NormalizationStrategy == null)
        {             
            badResults = new QuantificationResults
            {
                Summary = "QuantificationEngine Error: Normalization strategy is null."
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

        // 2) Roll up to peptides
        var peptideMatrix = RollUpStrategy.RollUpSpectralMatches(ExperimentalDesign, SpectralMatches, ModifiedBioPolymers);

        // 3) Normalize peptide matrix 
        var peptideMatrixNorm = NormalizationStrategy.NormalizePeptideIntensities(peptideMatrix, ModifiedBioPolymers);

        // 4) Roll up to protein groups
        var proteinMatrix = RollUpStrategy.RollUpPeptides(peptideMatrixNorm, BioPolymerGroups);

        // 5) Normalize protein group matrix
        var proteinMatrixNorm = NormalizationStrategy.NormalizeProteinIntensities(proteinMatrix, BioPolymerGroups);

        // 6) Write results (If enabled)
        if(Parameters.WritePeptideInformation)
            QuantificationWriter.WritePeptideMatrix(peptideMatrixNorm, Parameters.OutputDirectory);
        if(Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrixNorm, Parameters.OutputDirectory);

        return new QuantificationResults
        {
            Summary = "Quantification completed successfully."
        };
    }

    public QuantMatrix<IBioPolymerWithSetMods> Pivot(List<ISpectralMatch> spectralMatches, IExperimentalDesign experimentalDesign)
    {
        var sampleInfos = GetOrderedSampleInfos(experimentalDesign, spectralMatches, out var filePathToArrayPositionDict);
    }

    public static List<ISampleInfo> GetOrderedSampleInfos(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches)
    {
        var fileNames = spectralMatches
            .Select(sm => Path.GetFileName(sm.FullFilePath))
            .Distinct()
            .Order()
            .ToList();

        var samples = fileNames
            .SelectMany(fn => experimentalDesign.FileNameSampleInfoDictionary[fn])
            .ToList();

        return samples;
    }

    
}
