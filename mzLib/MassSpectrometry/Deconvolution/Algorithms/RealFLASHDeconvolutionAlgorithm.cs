using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Performs deconvolution via the official FLASHDeconv executable (OpenMS 3.0).
    ///
    /// Two usage modes:
    ///   1. Per-scan: Deconvolute(MzSpectrum, MzRange) — writes a temp single-scan mzML.
    ///      No multi-scan feature tracing; results may be sparse.
    ///   2. Whole-file: DeconvoluteFile(mzmlPath, params) — passes the full mzML directly.
    ///      Enables multi-scan feature tracing; strongly preferred.
    ///
    /// IMPORTANT: -out_spec requires one path per MS level (space-separated).
    ///   For MS1+MS2 files: -out_spec ms1.tsv ms2.tsv
    ///   For MS1-only files: -out_spec ms1.tsv
    ///   Passing only one path to a mixed file causes FLASHDeconv to write MS2 output
    ///   to that path and nothing for MS1, resulting in zero envelopes.
    ///   We always pass two paths and read the MS1 one.
    ///
    /// mzML requirements for OpenMS 3.0 pre-release:
    ///   - No UTF-8 BOM
    ///   - softwareList, instrumentConfigurationList (id=IC1), dataProcessingList (id=dp)
    ///   - run element:      defaultInstrumentConfigurationRef="IC1"
    ///   - spectrum element: defaultInstrumentConfigurationRef="IC1"
    ///
    /// TSV columns verified from real output:
    ///   ScanNum MonoisotopicMass SumIntensity RepresentativeCharge
    ///   RepresentativeMzStart RepresentativeMzEnd IsotopeCosine MassSNR Qscore
    /// </summary>
    internal class RealFLASHDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        /// <summary>
        /// Delegate signature for invoking FLASHDeconv. Production calls this
        /// with the real-Process implementation (<see cref="RunFLASHDeconvDefault"/>);
        /// tests inject a stub that writes canned TSV files to the requested output
        /// paths so the surrounding orchestration (mzML writing, range filtering,
        /// TSV parsing) can be unit-tested without launching a subprocess.
        /// </summary>
        internal delegate void FLASHDeconvRunner(
            string exePath,
            string inputMzml,
            string outFeatureTsv,
            string outMs1Tsv,
            string outMs2Tsv,
            RealFLASHDeconvolutionParameters parameters);

        /// <summary>
        /// Hardcoded install paths probed when no explicit FLASHDeconv path is supplied
        /// and PATH search misses. Exposed internally so tests can inject a known list
        /// (or an empty list) to exercise the resolver's branches deterministically.
        /// </summary>
        internal static readonly IReadOnlyList<string> DefaultWellKnownPaths = new[]
        {
            @"C:\Program Files\OpenMS-3.0.0-pre-HEAD-2023-06-17\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.5.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.4.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.3.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.0.0\bin\FLASHDeconv.exe",
            "/usr/bin/FLASHDeconv",
            "/usr/local/bin/FLASHDeconv",
            "/opt/openms/bin/FLASHDeconv",
        };

        private readonly FLASHDeconvRunner _runner;

        internal RealFLASHDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
            : this(deconParameters, runner: null) { }

        /// <summary>
        /// Test-friendly constructor that accepts an injected runner. The default
        /// (null) runner is the real Process-based implementation. Tests pass a
        /// stub that writes canned TSV files to the requested output paths.
        /// </summary>
        internal RealFLASHDeconvolutionAlgorithm(
            DeconvolutionParameters deconParameters,
            FLASHDeconvRunner? runner)
            : base(deconParameters)
        {
            _runner = runner ?? RunFLASHDeconvDefault;
        }

        // ── Per-scan entry point ──────────────────────────────────────────────

        /// <summary>
        /// Per-scan deconvolution. Range filtering uses the midpoint of FLASHDeconv's
        /// representative m/z range; envelopes whose isotope tail extends past
        /// spectrum.Range.Maximum may be dropped at the edge. For full-fidelity
        /// results use <see cref="DeconvoluteFile"/>.
        /// </summary>
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var p = DeconvolutionParameters as RealFLASHDeconvolutionParameters
                ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            range ??= spectrum.Range;
            if (spectrum.Size == 0) return Enumerable.Empty<IsotopicEnvelope>();

            string exePath = ResolveExePath(p.FLASHDeconvExePath);
            string guid = Guid.NewGuid().ToString("N");
            string tmpIn = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_in.mzML");
            string tmpFeat = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_feat.tsv");
            string tmpMs1 = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_ms1.tsv");
            string tmpMs2 = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_ms2.tsv");

            try
            {
                WriteSingleScanMzml(spectrum, tmpIn, p.Polarity);
                _runner(exePath, tmpIn, tmpFeat, tmpMs1, tmpMs2, p);
                if (!File.Exists(tmpMs1)) return Enumerable.Empty<IsotopicEnvelope>();
                return ParseSpecTsvByScan(tmpMs1, p)
                    .Values.SelectMany(x => x)
                    .Where(e => range.Contains(e.Peaks[0].mz));
            }
            finally
            {
                foreach (string f in Directory.GetFiles(p.WorkingDirectory, $"mzlib_fd_{guid}*"))
                    TryDelete(f);
            }
        }

        // ── Whole-file entry point ────────────────────────────────────────────

        /// <summary>
        /// Runs FLASHDeconv on a complete mzML file and returns all MS1 deconvolved
        /// envelopes grouped by 1-based scan number.
        /// This is the preferred mode — multi-scan feature tracing produces far better results.
        /// </summary>
        internal static Dictionary<int, List<IsotopicEnvelope>> DeconvoluteFile(
            string mzmlPath, RealFLASHDeconvolutionParameters p)
            => DeconvoluteFile(mzmlPath, p, runner: null);

        /// <summary>
        /// Test-friendly overload that accepts an injected <see cref="FLASHDeconvRunner"/>.
        /// The default (null) runner is the real Process-based implementation.
        /// </summary>
        internal static Dictionary<int, List<IsotopicEnvelope>> DeconvoluteFile(
            string mzmlPath,
            RealFLASHDeconvolutionParameters p,
            FLASHDeconvRunner? runner)
        {
            if (!File.Exists(mzmlPath))
                throw new FileNotFoundException($"mzML not found: {mzmlPath}", mzmlPath);

            FLASHDeconvRunner effectiveRunner = runner ?? RunFLASHDeconvDefault;
            string exePath = ResolveExePath(p.FLASHDeconvExePath);
            string guid = Guid.NewGuid().ToString("N");
            string tmpFeat = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_feat.tsv");
            string tmpMs1 = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_ms1.tsv");
            string tmpMs2 = Path.Combine(p.WorkingDirectory, $"mzlib_fd_{guid}_ms2.tsv");

            try
            {
                effectiveRunner(exePath, mzmlPath, tmpFeat, tmpMs1, tmpMs2, p);
                if (!File.Exists(tmpMs1)) return new Dictionary<int, List<IsotopicEnvelope>>();
                return ParseSpecTsvByScan(tmpMs1, p);
            }
            finally
            {
                foreach (string f in Directory.GetFiles(p.WorkingDirectory, $"mzlib_fd_{guid}*"))
                    TryDelete(f);
            }
        }

        // ── Exe resolution ────────────────────────────────────────────────────

        private static string ResolveExePath(string? explicitPath)
            => ResolveExePath(
                explicitPath,
                DefaultWellKnownPaths,
                Environment.GetEnvironmentVariable("PATH"));

        /// <summary>
        /// Test-friendly overload that takes the well-known-paths list and PATH-env
        /// value as parameters instead of reading them from compiled-in defaults +
        /// the live process environment. Lets tests deterministically exercise the
        /// well-known-search, PATH-search, and not-found-anywhere branches without
        /// depending on what's installed on the runner.
        /// </summary>
        internal static string ResolveExePath(
            string? explicitPath,
            IEnumerable<string> wellKnownPaths,
            string? pathEnv)
        {
            if (!string.IsNullOrWhiteSpace(explicitPath))
            {
                if (File.Exists(explicitPath)) return explicitPath!;
                throw new FileNotFoundException($"FLASHDeconv not found at: {explicitPath}", explicitPath);
            }

            foreach (string c in wellKnownPaths)
                if (File.Exists(c)) return c;

            if (pathEnv != null)
            {
                string[] names = OperatingSystem.IsWindows()
                    ? new[] { "FLASHDeconv.exe" }
                    : new[] { "FLASHDeconv", "flashdeconv" };
                foreach (string dir in pathEnv.Split(Path.PathSeparator))
                    foreach (string name in names)
                    {
                        string full = Path.Combine(dir, name);
                        if (File.Exists(full)) return full;
                    }
            }

            throw new FileNotFoundException(
                "FLASHDeconv not found. Set RealFLASHDeconvolutionParameters.FLASHDeconvExePath.");
        }

        // ── Process invocation ────────────────────────────────────────────────

        /// <summary>
        /// Default <see cref="FLASHDeconvRunner"/> implementation: spawns a real
        /// Process, drains its stdout/stderr asynchronously, enforces the configured
        /// timeout, and throws on non-zero exit.
        ///
        /// Excluded from code-coverage stats because every line is pure subprocess
        /// plumbing -- testing it requires either the real FLASHDeconv exe (covered
        /// by the RealFLASH_RealExe_* integration tests when OpenMS is installed)
        /// or a Process abstraction that's out of scope here. The orchestration
        /// around this method is testable via <see cref="FLASHDeconvRunner"/>
        /// injection on the constructor and on <see cref="DeconvoluteFile"/>.
        /// </summary>
        [ExcludeFromCodeCoverage]
        private static void RunFLASHDeconvDefault(
            string exePath, string inputMzml,
            string outFeatureTsv, string outMs1Tsv, string outMs2Tsv,
            RealFLASHDeconvolutionParameters p)
        {
            int minCharge = p.Polarity == Polarity.Negative
                ? -Math.Abs(p.MaxAssumedChargeState) : Math.Abs(p.MinAssumedChargeState);
            int maxCharge = p.Polarity == Polarity.Negative
                ? -Math.Abs(p.MinAssumedChargeState) : Math.Abs(p.MaxAssumedChargeState);

            var args = new StringBuilder();
            args.Append($"-in \"{inputMzml}\" ");
            args.Append($"-out \"{outFeatureTsv}\" ");
            // Two paths: MS1 first, MS2 second. Required for mixed-level files.
            args.Append($"-out_spec \"{outMs1Tsv}\" \"{outMs2Tsv}\" ");
            args.Append($"-Algorithm:tol {p.TolerancePpm.ToString(CultureInfo.InvariantCulture)} ");
            args.Append($"-Algorithm:min_charge {minCharge} ");
            args.Append($"-Algorithm:max_charge {maxCharge} ");
            args.Append($"-Algorithm:min_mass {p.MinMass.ToString(CultureInfo.InvariantCulture)} ");
            args.Append($"-Algorithm:max_mass {p.MaxMass.ToString(CultureInfo.InvariantCulture)} ");
            args.Append($"-Algorithm:min_isotope_cosine {p.MinIsotopeCosine.ToString(CultureInfo.InvariantCulture)} ");
            args.Append("-write_detail 0 -threads 1");

            using var proc = new Process();
            proc.StartInfo = new ProcessStartInfo
            {
                FileName = exePath,
                Arguments = args.ToString(),
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                WorkingDirectory = Path.GetDirectoryName(exePath) ?? ""
            };

            proc.Start();
            // Drain both streams asynchronously: a synchronous ReadToEnd before
            // WaitForExit defeats the timeout (it blocks until the child exits)
            // and risks a stdout/stderr pipe-buffer deadlock on chatty children.
            var stdoutTask = Task.Run(() => proc.StandardOutput.ReadToEnd());
            var stderrTask = Task.Run(() => proc.StandardError.ReadToEnd());

            bool done = proc.WaitForExit(p.ProcessTimeoutSeconds * 1_000);
            if (!done)
            {
                try { proc.Kill(); } catch { }
                throw new TimeoutException($"FLASHDeconv timed out after {p.ProcessTimeoutSeconds}s.");
            }

            string stdout = stdoutTask.Result;
            string stderr = stderrTask.Result;
            if (proc.ExitCode != 0)
                throw new MzLibException(
                    $"FLASHDeconv exited with code {proc.ExitCode}.\n" +
                    $"Args: {args}\nSTDOUT: {stdout}\nSTDERR: {stderr}");
        }

        // ── TSV parser ────────────────────────────────────────────────────────

        internal static Dictionary<int, List<IsotopicEnvelope>> ParseSpecTsvByScan(
            string tsvPath, RealFLASHDeconvolutionParameters p)
        {
            var results = new Dictionary<int, List<IsotopicEnvelope>>();
            string[] lines = File.ReadAllLines(tsvPath);
            if (lines.Length < 2) return results;

            // Trim trailing empty columns caused by a trailing tab in the header row
            string[] header = lines[0].Split('\t')
                .Where(h => !string.IsNullOrEmpty(h)).ToArray();

            int idxScan = Col(header, "ScanNum");
            int idxMono = Col(header, "MonoisotopicMass", "MonoMass");
            int idxZ = Col(header, "RepresentativeCharge", "RepCharge", "Charge");
            int idxInt = Col(header, "SumIntensity", "Intensity");
            int idxMzStart = Col(header, "RepresentativeMzStart", "RepMz", "Mz");
            int idxMzEnd = ColOpt(header, "RepresentativeMzEnd");
            int idxCosine = Col(header, "IsotopeCosine", "MinIsotopeCosine");
            int idxQscore = ColOpt(header, "Qscore", "Score");
            int idxSNR = ColOpt(header, "MassSNR", "SNR");

            for (int i = 1; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                if (string.IsNullOrEmpty(line)) continue;
                string[] c = line.Split('\t');
                if (c.Length < header.Length) continue;
                if (!I(c, idxScan, out int scanNum)) continue;
                if (!D(c, idxMono, out double mono)) continue;
                if (!I(c, idxZ, out int z)) continue;
                if (!D(c, idxInt, out double intensity)) continue;
                if (!D(c, idxMzStart, out double mzStart)) continue;
                double repMz = mzStart;
                if (idxMzEnd >= 0 && D(c, idxMzEnd, out double mzEnd))
                    repMz = (mzStart + mzEnd) / 2.0;
                D(c, idxCosine, out double cosine);
                if (!D(c, idxQscore, out double score))
                    if (!D(c, idxSNR, out score))
                        score = cosine;
                int signedZ = p.Polarity == Polarity.Negative ? -Math.Abs(z) : Math.Abs(z);
                if (mono < p.MinMass || mono > p.MaxMass) continue;
                var env = new IsotopicEnvelope(
                    id: i - 1,
                    peaks: new List<(double mz, double intensity)> { (repMz, intensity) },
                    monoisotopicmass: mono,
                    chargestate: signedZ,
                    intensity: intensity,
                    score: score);
                if (!results.TryGetValue(scanNum, out var list))
                    results[scanNum] = list = new List<IsotopicEnvelope>();
                list.Add(env);
            }
            return results;
        }

        // ── mzML writer ───────────────────────────────────────────────────────

        internal static void WriteSingleScanMzml(MzSpectrum spectrum, string path, Polarity polarity)
        {
            int n = spectrum.XArray.Length;
            double lowMz = n > 0 ? spectrum.XArray[0] : 0;
            double hiMz = n > 0 ? spectrum.XArray[n - 1] : 0;
            double tic = n > 0 ? spectrum.YArray.Sum() : 0;
            double bpInt = n > 0 ? spectrum.YArray.Max() : 0;
            double bpMz = n > 0 ? spectrum.XArray[Array.IndexOf(spectrum.YArray, bpInt)] : 0;

            string polTerm = polarity == Polarity.Negative
                ? @"<cvParam cvRef=""MS"" accession=""MS:1000129"" name=""negative scan""/>"
                : @"<cvParam cvRef=""MS"" accession=""MS:1000130"" name=""positive scan""/>";

            string mzB64 = B64(spectrum.XArray);
            string intB64 = B64(spectrum.YArray);

            var sb = new StringBuilder();
            sb.AppendLine(@"<?xml version=""1.0"" encoding=""utf-8""?>");
            sb.AppendLine(@"<mzML xmlns=""http://psi.hupo.org/ms/mzml"" xmlns:xsi=""http://www.w3.org/2001/XMLSchema-instance"" xsi:schemaLocation=""http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"">");
            sb.AppendLine(@"  <cvList count=""2"">");
            sb.AppendLine(@"    <cv id=""MS"" fullName=""Proteomics Standards Initiative Mass Spectrometry Ontology"" version=""4.1.30"" URI=""https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo""/>");
            sb.AppendLine(@"    <cv id=""UO"" fullName=""Unit Ontology"" URI=""https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo""/>");
            sb.AppendLine(@"  </cvList>");
            sb.AppendLine(@"  <fileDescription><fileContent><cvParam cvRef=""MS"" accession=""MS:1000579"" name=""MS1 spectrum""/></fileContent></fileDescription>");
            sb.AppendLine(@"  <softwareList count=""1""><software id=""mzLib"" version=""1.0""><cvParam cvRef=""MS"" accession=""MS:1000799"" name=""custom unreleased software tool"" value=""mzLib""/></software></softwareList>");
            sb.AppendLine(@"  <instrumentConfigurationList count=""1""><instrumentConfiguration id=""IC1""><cvParam cvRef=""MS"" accession=""MS:1000031"" name=""instrument model""/></instrumentConfiguration></instrumentConfigurationList>");
            sb.AppendLine(@"  <dataProcessingList count=""1""><dataProcessing id=""dp""><processingMethod order=""0"" softwareRef=""mzLib""><cvParam cvRef=""MS"" accession=""MS:1000544"" name=""Conversion to mzML""/></processingMethod></dataProcessing></dataProcessingList>");
            sb.AppendLine(@"  <run defaultInstrumentConfigurationRef=""IC1"">");
            sb.AppendLine(@"    <spectrumList count=""1"" defaultDataProcessingRef=""dp"">");
            sb.AppendLine($@"      <spectrum index=""0"" id=""scan=1"" defaultArrayLength=""{n}"" defaultInstrumentConfigurationRef=""IC1"">");
            sb.AppendLine(@"        <cvParam cvRef=""MS"" accession=""MS:1000511"" name=""ms level"" value=""1""/>");
            sb.AppendLine($"        {polTerm}");
            sb.AppendLine(@"        <cvParam cvRef=""MS"" accession=""MS:1000127"" name=""centroid spectrum""/>");
            sb.AppendLine(CV("MS:1000528", "lowest observed m/z", lowMz, "MS:1000040", "m/z"));
            sb.AppendLine(CV("MS:1000527", "highest observed m/z", hiMz, "MS:1000040", "m/z"));
            sb.AppendLine(CV("MS:1000285", "total ion current", tic));
            sb.AppendLine(CV("MS:1000504", "base peak m/z", bpMz, "MS:1000040", "m/z"));
            sb.AppendLine(CV("MS:1000505", "base peak intensity", bpInt));
            sb.AppendLine(@"        <scanList count=""1""><cvParam cvRef=""MS"" accession=""MS:1000795"" name=""no combination""/><scan><cvParam cvRef=""MS"" accession=""MS:1000016"" name=""scan start time"" value=""0"" unitCvRef=""UO"" unitAccession=""UO:0000031"" unitName=""minute""/></scan></scanList>");
            sb.AppendLine(@"        <binaryDataArrayList count=""2"">");
            sb.AppendLine($@"          <binaryDataArray encodedLength=""{mzB64.Length}""><cvParam cvRef=""MS"" accession=""MS:1000514"" name=""m/z array"" unitCvRef=""MS"" unitAccession=""MS:1000040"" unitName=""m/z""/><cvParam cvRef=""MS"" accession=""MS:1000576"" name=""no compression""/><cvParam cvRef=""MS"" accession=""MS:1000523"" name=""64-bit float""/><binary>{mzB64}</binary></binaryDataArray>");
            sb.AppendLine($@"          <binaryDataArray encodedLength=""{intB64.Length}""><cvParam cvRef=""MS"" accession=""MS:1000515"" name=""intensity array"" unitCvRef=""MS"" unitAccession=""MS:1000131"" unitName=""number of detector counts""/><cvParam cvRef=""MS"" accession=""MS:1000576"" name=""no compression""/><cvParam cvRef=""MS"" accession=""MS:1000523"" name=""64-bit float""/><binary>{intB64}</binary></binaryDataArray>");
            sb.AppendLine(@"        </binaryDataArrayList>");
            sb.AppendLine(@"      </spectrum>");
            sb.AppendLine(@"    </spectrumList>");
            sb.AppendLine(@"  </run>");
            sb.Append(@"</mzML>");

            File.WriteAllText(path, sb.ToString(),
                new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        }

        private static string CV(string acc, string name, double val, string? unitAcc = null, string? unitName = null)
        {
            string v = val.ToString(CultureInfo.InvariantCulture);
            string unit = unitAcc != null ? $@" unitCvRef=""MS"" unitAccession=""{unitAcc}"" unitName=""{unitName}""" : "";
            return $@"        <cvParam cvRef=""MS"" accession=""{acc}"" name=""{name}"" value=""{v}""{unit}/>";
        }

        private static string B64(double[] values)
        {
            var bytes = new byte[values.Length * 8];
            Buffer.BlockCopy(values, 0, bytes, 0, bytes.Length);
            return Convert.ToBase64String(bytes);
        }

        private static int Col(string[] header, params string[] names)
        {
            foreach (string name in names)
            {
                int idx = Array.FindIndex(header, h => string.Equals(h.Trim(), name, StringComparison.OrdinalIgnoreCase));
                if (idx >= 0) return idx;
            }
            throw new MzLibException(
                $"RealFLASHDeconvolution: required column not found. Tried: [{string.Join(", ", names)}]. Header: [{string.Join(", ", header)}]");
        }

        private static int ColOpt(string[] header, params string[] names)
        {
            foreach (string name in names)
            {
                int idx = Array.FindIndex(header, h => string.Equals(h.Trim(), name, StringComparison.OrdinalIgnoreCase));
                if (idx >= 0) return idx;
            }
            return -1;
        }

        private static bool D(string[] cols, int idx, out double v)
        {
            v = 0;
            if (idx < 0 || idx >= cols.Length) return false;
            return double.TryParse(cols[idx].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        private static bool I(string[] cols, int idx, out int v)
        {
            v = 0;
            if (idx < 0 || idx >= cols.Length) return false;
            return int.TryParse(cols[idx].Trim(), out v);
        }

        private static void TryDelete(string path)
        {
            try { if (File.Exists(path)) File.Delete(path); } catch { }
        }
    }
}