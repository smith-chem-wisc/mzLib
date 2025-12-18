using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;
using Quantification.Strategies;
using System;
using MassSpectrometry;
using MzLibUtil;

namespace Quantification;

/// <summary>
/// Orchestrates: (raw snapshot) -> write raw -> normalize -> rollups -> optional writes.
/// Does not depend on MetaMorpheus; MetaMorpheus supplies the objects.
/// </summary>
public sealed class QuantificationEngine
{
    public QuantificationParameters Parameters { get; init; }
    public NormalizationStrategyType NormalizationStrategyType { get; init; }
    public RollUpStrategyType RollUpStrategy { get; init; }
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
        NormalizationStrategyType = parameters.NormalizationStrategy;
        RollUpStrategy = parameters.RollUpStrategy;
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
        try { GetNormalizationStrategy(NormalizationStrategyType); }
        catch (NotImplementedException e)
        {
            badResults = new QuantificationResults
            {
                Summary = $"QuantificationEngine Error: {e.Message}"
            };
            return false;
        }
        try { GetRollUpStrategy(RollUpStrategy); }
        catch (NotImplementedException e)
        {
            badResults = new QuantificationResults
            {
                Summary = $"QuantificationEngine Error: {e.Message}"
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

        // 2) Create our normalization and roll-up strategies
        var normalization = GetNormalizationStrategy(NormalizationStrategyType);
        var rollup = GetRollUpStrategy(RollUpStrategy);

        // 3) Roll up to peptides
        var peptideMatrixRaw = rollup.RollUpSpectralMatches(ExperimentalDesign, SpectralMatches, ModifiedBioPolymers);

        // 4) Normalize peptide matrix 
        var peptideMatrixNorm = normalization.NormalizePeptideIntensities(peptideMatrixRaw, ModifiedBioPolymers);

        // 5) Roll up to protein groups
        var proteinMatrixRaw = rollup.RollUpPeptides(peptideMatrixNorm, BioPolymerGroups);

        // 6) Normalize protein group matrix
        var proteinMatrixNorm = normalization.NormalizeProteinIntensities(proteinMatrixRaw, BioPolymerGroups);

        // 7) Write results (If enabled)
        if(Parameters.WritePeptideInformation)
            QuantificationWriter.WritePeptideMatrix(peptideMatrixNorm, Parameters.OutputDirectory);
        if(Parameters.WriteProteinInformation)
            QuantificationWriter.WriteProteinGroupMatrix(proteinMatrixNorm, Parameters.OutputDirectory);

        return new();
    }

    private static INormalizationStrategy GetNormalizationStrategy(NormalizationStrategyType strategyType)
    {
        return strategyType switch
        {
            NormalizationStrategyType.None => new NoNormalization(),
            //NormalizationStrategyType.Median => new NormalizationStrategies.MedianNormalization(),
            //NormalizationStrategyType.Quantile => new NormalizationStrategies.QuantileNormalization(),
            //NormalizationStrategyType.VarianceStabilization => new NormalizationStrategies.VarianceStabilizationNormalization(),
            _ => throw new NotImplementedException($"Normalization strategy '{strategyType}' is not implemented."),
        };
    }

    private static IRollUpStrategy GetRollUpStrategy(RollUpStrategyType strategyType)
    {
        return strategyType switch
        {
            RollUpStrategyType.Sum => new SumRollUp(),
            //RollUpStrategyType.Median => new MedianRollUpStrategy(),
            //RollUpStrategyType.WeightedAverage => new WeightedAverageRollUpStrategy(),
            _ => throw new NotImplementedException($"Roll-up strategy '{strategyType}' is not implemented."),
        };
    }
}
