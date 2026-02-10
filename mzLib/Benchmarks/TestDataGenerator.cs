using Chromatography.RetentionTimePrediction;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace Benchmarks;

/// <summary>
/// Generates test peptides for benchmarking retention time predictors.
/// </summary>
public static class TestDataGenerator
{
    private static readonly string[] CommonSequences =
    [
        "PEPTIDE",
        "ALAVDGAGKPGAEE",
        "YHPDKNPSEEAAEK",
        "GVTVDKMTELR",
        "VGIVPGEVIAPGMR",
        "SEQUENCEWITHCYSTEINE",
        "ANOTHERTESTPEPTIDE",
        "VERYLONGPEPTIDESEQUENCEFORBENCHMARKING",
        "SHORTSEQ",
        "MEDIUMPEPTIDESEQ"
    ];

    /// <summary>
    /// Generates a list of unmodified peptides for benchmarking.
    /// </summary>
    public static List<IRetentionPredictable> GenerateUnmodifiedPeptides(int count)
    {
        var peptides = new List<IRetentionPredictable>(count);
        var emptyMods = new Dictionary<string, Modification>();

        for (int i = 0; i < count; i++)
        {
            var sequence = CommonSequences[i % CommonSequences.Length];
            peptides.Add(new PeptideWithSetModifications(sequence, emptyMods));
        }

        return peptides;
    }

    /// <summary>
    /// Generates peptides with varying sequence lengths.
    /// </summary>
    public static List<IRetentionPredictable> GeneratePeptidesWithVaryingLengths(int minLength, int maxLength, int count)
    {
        var peptides = new List<IRetentionPredictable>(count);
        var emptyMods = new Dictionary<string, Modification>();
        var random = new Random(42); // Fixed seed for reproducibility
        var aminoAcids = "ACDEFGHIKLMNPQRSTVWY";

        for (int i = 0; i < count; i++)
        {
            int length = minLength + (i % (maxLength - minLength + 1));
            var sequence = new string(Enumerable.Range(0, length)
                .Select(_ => aminoAcids[random.Next(aminoAcids.Length)])
                .ToArray());
            peptides.Add(new PeptideWithSetModifications(sequence, emptyMods));
        }

        return peptides;
    }
}