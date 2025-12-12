using Proteomics.ProteolyticDigestion;
namespace Proteomics.RetentionTimePrediction;

// This extension method and class ensure backwards compatibility with the old structure prior to the existence of the Chromatography namespace
public static class ChromatographyExtensions
{
    public static double ScoreSequence(this Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3 predictor, PeptideWithSetModifications peptide) => predictor.ScoreSequence(peptide.BaseSequence);
}

public class SSRCalc3: Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3
{
    public SSRCalc3(string name, Column column) : base(name, column) { }
}