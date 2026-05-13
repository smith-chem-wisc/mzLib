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
            Results = csv.GetRecords<Ms1Feature>().ToList();

            Software = Results.All(p => p.IntensityApex == null) ? Software.FLASHDeconv : Software.TopFD;
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
        /// end/apex, charge bounds, apex intensity) survives. The schema is
        /// FLASHDeconv-style (newer-TopFD column aliases are aliases, not
        /// the canonical write target).
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
        /// <param name="software">Software label for the produced file.
        /// Defaults to <see cref="Software.FLASHDeconv"/> because that's
        /// the canonical write schema; pass <see cref="Software.Unspecified"/>
        /// to leave it blank.</param>
        public static Ms1FeatureFile FromMassFeatures(
            IEnumerable<MassFeature> features,
            int sampleId = 0,
            int fractionId = 0,
            Software software = Software.FLASHDeconv)
        {
            var file = new Ms1FeatureFile { Software = software };
            int sequentialId = 0;
            foreach (var f in features)
            {
                file.Results.Add(f.ToMs1Feature(sequentialId++, sampleId, fractionId));
            }
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
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
