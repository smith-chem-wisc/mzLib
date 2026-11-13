using Chromatography.RetentionTimePrediction;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Readers.DeNovoMapping;

/// <summary>
/// Settings for mapping de novo (Casanovo) predictions onto protein databases.
/// </summary>
public class DeNovoMappingParameters
{
    /// <summary>
    /// How database proteins are digested before their peptides are compared to the predictions.
    /// Defaults to trypsin, two missed cleavages, peptides of 7 residues or longer.
    /// </summary>
    public DigestionParams DigestionParams { get; init; } = new(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7);

    /// <summary>
    /// Decoys generated when a database is read from disk. Ignored when proteins are passed in directly;
    /// the caller is then expected to include its own decoys.
    /// </summary>
    public DecoyType DecoyType { get; init; } = DecoyType.Reverse;

    /// <summary>
    /// Treat isoleucine and leucine as the same residue. Both the predicted and the database peptides are
    /// converted to leucine before comparison, so the result does not depend on whether Casanovo was run
    /// with replace_isoleucine_with_leucine.
    /// </summary>
    public bool ReplaceIWithL { get; init; } = true;

    /// <summary>
    /// Fill in observed retention times from the spectra files (when the mzTab lacks them) and predict
    /// retention times with Chronologer, so each database gets a retention-time error distribution.
    /// </summary>
    public bool MapRetentionTimes { get; init; } = true;

    /// <summary>
    /// Spectra files that Casanovo was run on. They are matched to each mzTab ms_run by file name without
    /// extension, so the files need not sit where the mzTab says they were. Only read when an mzTab
    /// has PSMs without a retention time.
    /// </summary>
    public IReadOnlyList<string> SpectraFilePaths { get; init; } = [];

    /// <summary>
    /// Predictor for the expected retention time. Defaults to Chronologer, which is what MetaMorpheus uses for
    /// its own retention-time metrics. Earlier versions of this mapper (MetaCurator) used SSRCalc3.
    /// </summary>
    public RetentionTimePredictor? RetentionTimePredictor { get; init; }

    public int MaxThreads { get; init; } = Environment.ProcessorCount;
}
