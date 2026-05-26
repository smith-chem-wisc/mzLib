using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MassSpectrometry.Deconvolution.FeatureTracing;

namespace Readers
{
    /// <summary>
    /// Whole-file driver that turns a spectra file into a FLASHDeconv-style
    /// <c>_ms1.feature</c> file: deconvolve every MS1 scan with the supplied
    /// <see cref="DeconvolutionParameters"/> (e.g. <c>MetaFlashDeconParameters</c>), trace
    /// cross-scan mass features with the supplied <see cref="IMassFeatureTracer"/>
    /// (MetaFlashDecon-native or consensus), and write the result via
    /// <see cref="Ms1FeatureFile.FromMassFeatures"/>. This is the entry point that makes
    /// MetaFlashDecon "analogous to #1069": whole-file deconvolution &#8594; <c>_ms1.feature</c>
    /// &#8594; MetaMorpheus FromFile search.
    ///
    /// <para>
    /// <b>Retention time is written in SECONDS</b> to match real OpenMS FLASHDeconv output
    /// (so file-vs-file comparison is direct; MetaMorpheus's FromFile reader auto-normalises
    /// either unit). mzLib <see cref="MsDataScan.RetentionTime"/> is in minutes, so RTs are
    /// scaled by <c>secondsPerRtUnit</c> (default 60). Pass 1.0 if the source RT is already
    /// in seconds.
    /// </para>
    /// <para>
    /// <b>Schema:</b> the written file is the TopFD-flavoured <c>_ms1.feature</c> — the 11
    /// canonical FLASHDeconv columns (Sample_ID, ID, Mass, Intensity, Time_begin, Time_end,
    /// Time_apex, Minimum/Maximum_charge_state, Minimum/Maximum_fraction_id) PLUS an
    /// <c>Apex_intensity</c> column. It is a superset of the bare FLASHDeconv schema, read as
    /// the same <c>_ms1.feature</c> type by <see cref="Ms1FeatureFile"/> and MetaMorpheus
    /// (which uses Apex_intensity for per-charge weighting).
    /// </para>
    /// </summary>
    public static class Ms1FeatureGenerator
    {
        /// <summary>mzML <see cref="MsDataScan.RetentionTime"/> is in minutes; ×60 writes seconds.</summary>
        public const double DefaultSecondsPerRtUnit = 60.0;

        /// <summary>
        /// Deconvolve every MS1 scan and trace cross-scan mass features. Returned features carry
        /// retention times in the source's native units (minutes for mzML); RT-to-seconds
        /// conversion happens at write time (<see cref="WriteFeatures"/>).
        /// </summary>
        public static List<MassFeature> GenerateFeatures(
            string spectraFilePath,
            DeconvolutionParameters deconParameters,
            IMassFeatureTracer tracer)
        {
            var dataFile = MsDataFileReader.GetDataFile(spectraFilePath).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToList();

            var perScanEnvelopes = ms1Scans
                .Select(s => (IReadOnlyList<IsotopicEnvelope>)
                    Deconvoluter.Deconvolute(s.MassSpectrum, deconParameters).ToList())
                .ToList();

            return tracer.TraceFeatures(ms1Scans, perScanEnvelopes);
        }

        /// <summary>
        /// Write traced features to an <c>_ms1.feature</c> file, converting retention times to
        /// seconds. NOTE: this mutates the supplied features' envelope retention times in place
        /// (scaled by <paramref name="secondsPerRtUnit"/>, then re-finalised).
        /// </summary>
        public static void WriteFeatures(
            IReadOnlyList<MassFeature> features,
            string outputMs1FeaturePath,
            double secondsPerRtUnit = DefaultSecondsPerRtUnit)
        {
            ConvertRetentionTimeToSeconds(features, secondsPerRtUnit);
            Ms1FeatureFile.FromMassFeatures(features).WriteResults(outputMs1FeaturePath);
        }

        /// <summary>
        /// End-to-end: deconvolve + trace a spectra file and write its <c>_ms1.feature</c>
        /// (RT in seconds). Returns the features (retention times converted to seconds).
        /// </summary>
        public static List<MassFeature> GenerateAndWrite(
            string spectraFilePath,
            DeconvolutionParameters deconParameters,
            IMassFeatureTracer tracer,
            string outputMs1FeaturePath,
            double secondsPerRtUnit = DefaultSecondsPerRtUnit)
        {
            var features = GenerateFeatures(spectraFilePath, deconParameters, tracer);
            WriteFeatures(features, outputMs1FeaturePath, secondsPerRtUnit);
            return features;
        }

        private static void ConvertRetentionTimeToSeconds(IReadOnlyList<MassFeature> features, double scale)
        {
            if (scale == 1.0) return;
            foreach (var f in features)
            {
                foreach (var trace in f.Traces)
                    foreach (var e in trace.Envelopes)
                        e.RT *= scale;
                f.Finalise(); // recompute RTStart/RTEnd from the scaled envelope RTs
            }
        }
    }
}
