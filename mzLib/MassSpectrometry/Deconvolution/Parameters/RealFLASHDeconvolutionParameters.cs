#nullable enable
using System;
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// Parameters for <see cref="RealFLASHDeconvolutionAlgorithm"/>, which wraps
    /// the official FLASHDeconv executable from OpenMS.
    ///
    /// The algorithm does no path discovery at decon time. Register the
    /// executable once per process with <see cref="FlashDeconvExePathRegistry.Register"/>
    /// (optionally feeding it <see cref="FlashDeconvExePathRegistry.Resolve"/>,
    /// which checks well-known install locations and PATH) and leave
    /// <see cref="FLASHDeconvExePath"/> null; or set <see cref="FLASHDeconvExePath"/>
    /// to override the registered path for this parameters object.
    /// </summary>
    public class RealFLASHDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.RealFLASHDeconvolution;

        // ── Executable location ───────────────────────────────────────────────

        /// <summary>
        /// Optional per-parameters override of the FLASHDeconv path. When null the
        /// algorithm uses <see cref="FlashDeconvExePathRegistry.RegisteredPath"/>,
        /// and throws if nothing is registered. A path set here is validated once
        /// per process (cached by the registry), not on every decon call.
        /// </summary>
        public string? FLASHDeconvExePath { get; set; }

        /// <summary>
        /// Directory for temporary mzML / TSV files.  Defaults to system temp.
        /// All temp files are deleted after each call.
        /// </summary>
        public string WorkingDirectory { get; set; } = Path.GetTempPath();

        /// <summary>
        /// Seconds to wait before killing the FLASHDeconv process.
        /// Default: 300 s (5 minutes).
        /// </summary>
        public int ProcessTimeoutSeconds { get; set; } = 300;

        // ── Algorithm flags (map to FLASHDeconv CLI) ──────────────────────────

        /// <summary>ppm tolerance  →  -Algorithm:tol</summary>
        public double TolerancePpm { get; set; } = 10.0;

        /// <summary>Minimum neutral mass (Da)  →  -Algorithm:min_mass</summary>
        public double MinMass { get; set; } = 50.0;

        /// <summary>Maximum neutral mass (Da)  →  -Algorithm:max_mass</summary>
        public double MaxMass { get; set; } = 100_000.0;

        /// <summary>
        /// Minimum isotope cosine similarity  →  -Algorithm:min_isotope_cosine
        /// Default matches FLASHDeconv's own default (0.85).
        /// </summary>
        public double MinIsotopeCosine { get; set; } = 0.85;

        // ── Constructor ───────────────────────────────────────────────────────

        public RealFLASHDeconvolutionParameters(
            int minCharge = 1,
            int maxCharge = 60,
            double tolerancePpm = 10.0,
            double minMass = 50.0,
            double maxMass = 100_000.0,
            double minIsotopeCosine = 0.85,
            Polarity polarity = Polarity.Positive,
            string? flashDeconvExePath = null,
            string? workingDirectory = null,
            int processTimeoutSeconds = 300,
            AverageResidue? averageResidueModel = null)
            : base(minCharge, maxCharge, polarity, averageResidueModel)
        {
            TolerancePpm = tolerancePpm;
            MinMass = minMass;
            MaxMass = maxMass;
            MinIsotopeCosine = minIsotopeCosine;
            FLASHDeconvExePath = flashDeconvExePath;
            WorkingDirectory = workingDirectory ?? Path.GetTempPath();
            ProcessTimeoutSeconds = processTimeoutSeconds;
        }

        // Decoy deconvolution doesn't apply to the FLASHDeconv exe wrapper.
        public override DeconvolutionParameters? ToDecoyParameters() => null;

        #region IEquatable<RealFLASHDeconvolutionParameters>

        protected override bool EqualProperties(DeconvolutionParameters other)
        {
            var o = (RealFLASHDeconvolutionParameters)other;
            return TolerancePpm.Equals(o.TolerancePpm)
                && MinMass.Equals(o.MinMass)
                && MaxMass.Equals(o.MaxMass)
                && MinIsotopeCosine.Equals(o.MinIsotopeCosine)
                && string.Equals(FLASHDeconvExePath, o.FLASHDeconvExePath, StringComparison.Ordinal)
                && string.Equals(WorkingDirectory, o.WorkingDirectory, StringComparison.Ordinal)
                && ProcessTimeoutSeconds == o.ProcessTimeoutSeconds;
        }

        protected override void AddHashCodes(HashCode hash)
        {
            hash.Add(TolerancePpm);
            hash.Add(MinMass);
            hash.Add(MaxMass);
            hash.Add(MinIsotopeCosine);
            hash.Add(FLASHDeconvExePath, StringComparer.Ordinal);
            hash.Add(WorkingDirectory, StringComparer.Ordinal);
            hash.Add(ProcessTimeoutSeconds);
        }

        public override RealFLASHDeconvolutionParameters Clone()
        {
            return new RealFLASHDeconvolutionParameters(
                MinAssumedChargeState, MaxAssumedChargeState,
                TolerancePpm, MinMass, MaxMass, MinIsotopeCosine,
                Polarity, FLASHDeconvExePath, WorkingDirectory,
                ProcessTimeoutSeconds, AverageResidueModel)
            {
                UseGenericScore = UseGenericScore,
                ExpectedIsotopeSpacing = ExpectedIsotopeSpacing
            };
        }

        #endregion
    }
}