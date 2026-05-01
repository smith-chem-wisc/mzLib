#nullable enable
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// Parameters for <see cref="RealFLASHDeconvolutionAlgorithm"/>, which wraps
    /// the official FLASHDeconv executable from OpenMS.
    ///
    /// The executable is located via:
    ///   1. <see cref="FLASHDeconvExePath"/> (explicit, highest priority)
    ///   2. Well-known OpenMS install paths checked at runtime
    ///   3. System PATH
    ///
    /// In MetaMorpheus, set <c>GlobalSettings.FLASHDeconvExecutablePath</c> once
    /// (persisted in GlobalSettings.toml) and leave <see cref="FLASHDeconvExePath"/>
    /// null — the algorithm will use the global setting automatically.
    /// </summary>
    public class RealFLASHDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.RealFLASHDeconvolution;

        // ── Executable location ───────────────────────────────────────────────

        /// <summary>
        /// Full path to FLASHDeconv.exe.  Leave null to use the well-known-path
        /// search or GlobalSettings.FLASHDeconvExecutablePath.
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
    }
}