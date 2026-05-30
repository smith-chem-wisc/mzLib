using CsvHelper;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a ms1.feature deconvolution result
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class Ms1FeatureFile : ResultFile<Ms1Feature>, IResultFile, IMs1FeatureFile
    {
        public override SupportedFileType FileType => SupportedFileType.Ms1Feature;

        public sealed override Software Software { get; set; }

        public IEnumerable<ISingleChargeMs1Feature> GetMs1Features() => Results.SelectMany(r => r.GetSingleChargeFeatures());

        public Ms1FeatureFile(string filePath, Software deconSoftware = Software.Unspecified) : base(filePath,
            deconSoftware)
        {
            using (var sr = new StreamReader(filePath))
            {
                string firstLine = sr.ReadLine() ?? "";
                if (firstLine.Contains("\tApex_intensity\t") || firstLine.Contains("\tIntensity_Apex\t"))
                    Software = Software.TopFD;
                else
                    Software = Software.FLASHDeconv;
            }
        }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public Ms1FeatureFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), Ms1Feature.CsvConfiguration);
            var loaded = csv.GetRecords<Ms1Feature>().ToList();
            Results = loaded;

            // Use the just-loaded list, NOT the Results getter: for a header-only (zero-row)
            // file the getter sees an empty _results and re-invokes LoadResults, recursing
            // until the stack overflows.
            Software = loaded.All(p => p.IntensityApex == null) ? Software.FLASHDeconv : Software.TopFD;
        }

        /// <summary>
        /// Build an <see cref="Ms1FeatureFile"/> from a sequence of consensus-tracer
        /// <see cref="MassFeature"/> objects, pre-populated and ready to write
        /// via <see cref="WriteResults"/>. No file I/O happens in this factory
        /// -- the caller picks the destination path and invokes WriteResults
        /// separately.
        ///
        /// The output round-trips through <see cref="LoadResults"/>: every
        /// downstream-consumed field (Mass, Intensity, RetentionTime begin/
        /// end/apex, charge bounds, apex intensity) survives. Software is not
        /// stored as a column; it is re-derived on read from the presence of
        /// an Apex_intensity column, which the writer always emits -- so a
        /// reloaded factory file is detected as <see cref="Software.TopFD"/>.
        /// </summary>
        /// <param name="features">Finalised cross-charge features to emit.
        /// Caller must have invoked <see cref="MassFeature.Finalise"/> on
        /// each.</param>
        /// <param name="sampleId">Sample identifier written into every row.
        /// Default 0 fits the single-mzML case; multi-file producers can
        /// pass an external index.</param>
        /// <param name="fractionId">Fraction identifier written into both the
        /// Minimum_fraction_id and Maximum_fraction_id columns of every row.
        /// Default 0 fits a single fraction.</param>
        /// <param name="software">In-memory software label for the produced file.
        /// Defaults to <see cref="Software.TopFD"/> to match how the file is
        /// re-detected on reload (the writer emits an Apex_intensity column); the
        /// label is not persisted as a column. Pass <see cref="Software.Unspecified"/>
        /// to leave the in-memory label blank.</param>
        public static Ms1FeatureFile FromMassFeatures(
            IEnumerable<MassFeature> features,
            int sampleId = 0,
            int fractionId = 0,
            Software software = Software.TopFD)
        {
            // Build the record list up front and assign via the setter. The base
            // Results getter only lazy-loads when File.Exists(FilePath), so a
            // factory-built file (empty FilePath) returns these set records directly
            // -- single source of truth, no separate in-memory store needed.
            var records = new List<Ms1Feature>();
            int sequentialId = 0;
            foreach (var f in features)
            {
                // ONE row per feature with the aggregated Min/Max charge — faithful to the real
                // OpenMS FLASHDeconv _ms1.feature (FLASHDeconvFeatureFile writes a single row per
                // mass feature; VERIFIED 2026-05-30: in the ground-truth file every multi-row mass
                // differs in RETENTION TIME, never just charge — 0 charge-only splits). The older
                // ToMs1Features charge-gap splitting was unfaithful AND inflated intensity (it wrote
                // the full feature.SummedIntensity on every split row). Re-expanding Min..Max on
                // reload may fabricate an unobserved intermediate charge, but that matches FLASHDeconv
                // and the charge set is metadata, not used for mass/RT pairing.
                records.Add(f.ToMs1Feature(sequentialId++, sampleId, fractionId));
            }
            var file = new Ms1FeatureFile { Software = software };
            file.Results = records;
            return file;
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), Ms1Feature.CsvConfiguration);

            csv.WriteHeader<Ms1Feature>();
            // Results returns the in-memory factory records as-is (may be empty -> header
            // only), or lazy-loads from disk for a file-backed instance.
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
