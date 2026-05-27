using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// EXPLICIT differential tests: run the SAME fixed input through our C# port and through a
    /// standalone MSVC-compiled extract of the real OpenMS function, and compare. The C++ harness
    /// (E:\CodeReview\MetaFlashDecon\difftest\) writes its input to a shared file; this test reads
    /// that file so both sides see identical data. This is the anti-guessing check: it pins each
    /// ported function to the reference implementation's exact numeric output.
    /// </summary>
    [TestFixture]
    [Explicit("Differential test vs MSVC-compiled OpenMS extract; needs E:\\CodeReview\\MetaFlashDecon\\difftest input files.")]
    [ExcludeFromCodeCoverage]
    public class MetaFlashDeconDiffTest
    {
        private const string DiffDir = @"E:\CodeReview\MetaFlashDecon\difftest";

        [Test]
        public void NoisePeakPower_MatchesOpenMS()      => RunCase("noise_input.txt",     "noise_cpp_result.txt",     "noise_cs_result.txt");

        [Test]
        public void NoisePeakPower_Big_MatchesOpenMS()  => RunCase("noise_input_big.txt", "noise_cpp_result_big.txt", "noise_cs_result_big.txt");

        [Test] public void GetCosine1_MatchesOpenMS() => RunCosineCase(1);
        [Test] public void GetCosine2_MatchesOpenMS() => RunCosineCase(2);
        [Test] public void GetCosine3_MatchesOpenMS() => RunCosineCase(3);

        [Test]
        public void Recruit_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "recruit_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run recruit_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            double monoMass = double.Parse(h[0], CultureInfo.InvariantCulture);
            double isoDa = double.Parse(h[1], CultureInfo.InvariantCulture);
            int minCharge = int.Parse(h[2], CultureInfo.InvariantCulture);
            int maxCharge = int.Parse(h[3], CultureInfo.InvariantCulture);
            int minNegIso = int.Parse(h[4], CultureInfo.InvariantCulture);
            double tol = double.Parse(h[5], CultureInfo.InvariantCulture);
            int minIsotope = int.Parse(h[6], CultureInfo.InvariantCulture);
            int maxIsotope = int.Parse(h[7], CultureInfo.InvariantCulture);

            var mz = new List<double>();
            var inten = new List<double>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                mz.Add(double.Parse(t[0], CultureInfo.InvariantCulture));
                inten.Add(double.Parse(t[1], CultureInfo.InvariantCulture));
            }
            var spectrum = new MzSpectrum(mz.ToArray(), inten.ToArray(), false);

            var pg = new MetaFlashDeconPeakGroup
            {
                MinAbsCharge = minCharge,
                MaxAbsCharge = maxCharge,
                IsoDaDistance = isoDa,
                MinNegativeIsotopeIndex = minNegIso,
                IsPositive = true,
            };
            var noisy = pg.RecruitAllPeaksInSpectrum(spectrum, tol, monoMass, minIsotope, maxIsotope, Polarity.Positive);

            var recs = new List<(char kind, int charge, int iso, double mz, double intensity)>();
            foreach (var p in pg.SignalPeaks)      recs.Add(('S', p.AbsCharge, p.IsotopeIndex, p.Mz, p.Intensity));
            foreach (var p in pg.NegativeIsoPeaks) recs.Add(('G', p.AbsCharge, p.IsotopeIndex, p.Mz, p.Intensity));
            foreach (var p in noisy)               recs.Add(('N', p.AbsCharge, p.IsotopeIndex, p.Mz, p.Intensity));
            recs.Sort((a, b) =>
            {
                int k = a.kind.CompareTo(b.kind); if (k != 0) return k;
                if (a.charge != b.charge) return a.charge.CompareTo(b.charge);
                if (a.iso != b.iso) return a.iso.CompareTo(b.iso);
                return a.mz.CompareTo(b.mz);
            });
            var csLines = new List<string>();
            foreach (var r in recs)
                csLines.Add($"{r.kind} {r.charge} {r.iso} {r.mz.ToString("F6", CultureInfo.InvariantCulture)} {r.intensity.ToString("F6", CultureInfo.InvariantCulture)}");

            string csText = string.Join("\n", csLines);
            File.WriteAllText(Path.Combine(DiffDir, "recruit_cs_result.txt"), csText + "\n");
            TestContext.Progress.WriteLine($"CS recruit: {recs.Count} peaks");

            string cppResultPath = Path.Combine(DiffDir, "recruit_cpp_result.txt");
            if (File.Exists(cppResultPath))
            {
                var cppLines = new List<string>();
                foreach (var l in File.ReadAllLines(cppResultPath))
                    if (!string.IsNullOrWhiteSpace(l)) cppLines.Add(l.Trim());
                Assert.That(csLines.Count, Is.EqualTo(cppLines.Count),
                    $"peak count C# {csLines.Count} != OpenMS {cppLines.Count}");
                for (int i = 0; i < cppLines.Count; i++)
                    Assert.That(csLines[i], Is.EqualTo(cppLines[i]), $"recruited peak {i} differs");
            }
        }

        private static void RunCosineCase(int idx)
        {
            string inputPath = Path.Combine(DiffDir, $"cosine_input_{idx}.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run cosine_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int aStart = int.Parse(h[0], CultureInfo.InvariantCulture);
            int aEnd = int.Parse(h[1], CultureInfo.InvariantCulture);
            int offset = int.Parse(h[2], CultureInfo.InvariantCulture);
            int minIso = int.Parse(h[3], CultureInfo.InvariantCulture);
            int bSize = int.Parse(h[4], CultureInfo.InvariantCulture);
            double[] a = ParseRow(lines[1]);
            double[] b = ParseRow(lines[2]);

            double cs = MetaFlashDeconPeakGroup.GetCosine(a, aStart, aEnd, b, bSize, offset, minIso);
            TestContext.Progress.WriteLine($"CS  getCosine[{idx}] = {cs:F8}");

            string cppResultPath = Path.Combine(DiffDir, $"cosine_cpp_result_{idx}.txt");
            if (File.Exists(cppResultPath))
            {
                double cpp = double.Parse(File.ReadAllText(cppResultPath).Trim(), CultureInfo.InvariantCulture);
                Assert.That(cs, Is.EqualTo(cpp).Within(1e-5), $"C# ({cs}) != OpenMS ({cpp}) for getCosine[{idx}]");
            }
        }

        [Test]
        public void PerChargeInfo_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "percharge_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run perchargesnr_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minC = int.Parse(h[0], CultureInfo.InvariantCulture);
            int maxC = int.Parse(h[1], CultureInfo.InvariantCulture);
            double isoDa = double.Parse(h[2], CultureInfo.InvariantCulture);

            var pg = new MetaFlashDeconPeakGroup { MinAbsCharge = minC, MaxAbsCharge = maxC, IsoDaDistance = isoDa };
            var noisy = new List<MetaFlashDeconAlgorithm.LogMzPeak>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                double mz = double.Parse(t[1], CultureInfo.InvariantCulture);
                double inten = double.Parse(t[2], CultureInfo.InvariantCulture);
                int z = int.Parse(t[3], CultureInfo.InvariantCulture);
                var p = new MetaFlashDeconAlgorithm.LogMzPeak(mz, inten, 0.0, z, 0);
                if (t[0] == "S") pg.SignalPeaks.Add(p); else noisy.Add(p);
            }

            pg.UpdatePerChargeInformation(noisy);

            foreach (var l in File.ReadAllLines(Path.Combine(DiffDir, "percharge_cpp_result.txt")))
            {
                if (string.IsNullOrWhiteSpace(l)) continue;
                var t = l.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                int z = int.Parse(t[0], CultureInfo.InvariantCulture);
                AssertClose(pg.PerChargeInt[z], t[1], $"int z{z}");
                AssertClose(pg.PerChargeSumSignalSq[z], t[2], $"sumSq z{z}");
                AssertClose(pg.PerChargeNoisePwr[z], t[3], $"noise z{z}");
            }
        }

        [Test]
        public void Snr_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "snr_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run perchargesnr_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minC = int.Parse(h[0], CultureInfo.InvariantCulture);
            int maxC = int.Parse(h[1], CultureInfo.InvariantCulture);
            double globalCos = double.Parse(h[2], CultureInfo.InvariantCulture);

            int sz = 1 + maxC;
            var pg = new MetaFlashDeconPeakGroup
            {
                MinAbsCharge = minC, MaxAbsCharge = maxC, IsotopeCosineScore = globalCos,
                PerChargeInt = new double[sz], PerChargeSumSignalSq = new double[sz],
                PerChargeNoisePwr = new double[sz], PerChargeCos = new double[sz],
            };
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                int z = int.Parse(t[0], CultureInfo.InvariantCulture);
                pg.PerChargeInt[z] = double.Parse(t[1], CultureInfo.InvariantCulture);
                pg.PerChargeSumSignalSq[z] = double.Parse(t[2], CultureInfo.InvariantCulture);
                pg.PerChargeNoisePwr[z] = double.Parse(t[3], CultureInfo.InvariantCulture);
                pg.PerChargeCos[z] = double.Parse(t[4], CultureInfo.InvariantCulture);
            }

            pg.UpdateSNR();
            TestContext.Progress.WriteLine($"CS snr = {pg.Snr:F8}");

            foreach (var l in File.ReadAllLines(Path.Combine(DiffDir, "snr_cpp_result.txt")))
            {
                if (string.IsNullOrWhiteSpace(l)) continue;
                var t = l.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (t[0] == "snr") AssertClose(pg.Snr, t[1], "overall snr");
                else AssertClose(pg.GetChargeSnr(int.Parse(t[0], CultureInfo.InvariantCulture)), t[1], $"snr z{t[0]}");
            }
        }

        [Test]
        public void ChargeRange_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "chargerange_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run chain_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minC = int.Parse(h[0], CultureInfo.InvariantCulture);
            int maxC = int.Parse(h[1], CultureInfo.InvariantCulture);
            int sz = 1 + maxC;
            var pg = new MetaFlashDeconPeakGroup
            { MinAbsCharge = minC, MaxAbsCharge = maxC, PerChargeInt = new double[sz], PerChargeNoisePwr = new double[sz] };
            var noisy = new List<MetaFlashDeconAlgorithm.LogMzPeak>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (t[0] == "I") { int z = int.Parse(t[1], CultureInfo.InvariantCulture); pg.PerChargeInt[z] = double.Parse(t[2], CultureInfo.InvariantCulture); pg.PerChargeNoisePwr[z] = double.Parse(t[3], CultureInfo.InvariantCulture); }
                else { var p = new MetaFlashDeconAlgorithm.LogMzPeak(double.Parse(t[1], CultureInfo.InvariantCulture), double.Parse(t[2], CultureInfo.InvariantCulture), 0.0, int.Parse(t[3], CultureInfo.InvariantCulture), 0); if (t[0] == "S") pg.SignalPeaks.Add(p); else noisy.Add(p); }
            }
            pg.UpdateChargeRange(noisy);
            var cpp = File.ReadAllText(Path.Combine(DiffDir, "chargerange_cpp_result.txt")).Split(' ', StringSplitOptions.RemoveEmptyEntries);
            Assert.That(pg.MinAbsCharge, Is.EqualTo(int.Parse(cpp[0], CultureInfo.InvariantCulture)), "newMin");
            Assert.That(pg.MaxAbsCharge, Is.EqualTo(int.Parse(cpp[1], CultureInfo.InvariantCulture)), "newMax");
            Assert.That(pg.SignalPeaks.Count, Is.EqualTo(int.Parse(cpp[2], CultureInfo.InvariantCulture)), "filtered signal");
            Assert.That(noisy.Count, Is.EqualTo(int.Parse(cpp[3], CultureInfo.InvariantCulture)), "filtered noisy");
        }

        [Test]
        public void ChargeFit_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "chargefit_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run chain_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minC = int.Parse(h[0], CultureInfo.InvariantCulture);
            int maxC = int.Parse(h[1], CultureInfo.InvariantCulture);
            var pg = new MetaFlashDeconPeakGroup { MinAbsCharge = minC, MaxAbsCharge = maxC, PerChargeInt = new double[1 + maxC] };
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                pg.PerChargeInt[int.Parse(t[0], CultureInfo.InvariantCulture)] = double.Parse(t[1], CultureInfo.InvariantCulture);
            }
            pg.UpdateChargeFitScoreAndChargeIntensities();
            AssertClose(pg.ChargeScore, File.ReadAllText(Path.Combine(DiffDir, "chargefit_cpp_result.txt")).Trim(), "chargeScore");
        }

        [Test]
        public void MonoMass_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "monomass_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run chain_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minNegIso = int.Parse(h[0], CultureInfo.InvariantCulture);
            double isoDa = double.Parse(h[1], CultureInfo.InvariantCulture);
            var pg = new MetaFlashDeconPeakGroup { MinNegativeIsotopeIndex = minNegIso, IsoDaDistance = isoDa, IsPositive = true };
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                var p = new MetaFlashDeconAlgorithm.LogMzPeak(double.Parse(t[1], CultureInfo.InvariantCulture), double.Parse(t[2], CultureInfo.InvariantCulture), 0.0, int.Parse(t[3], CultureInfo.InvariantCulture), int.Parse(t[4], CultureInfo.InvariantCulture));
                if (t[0] == "S") pg.SignalPeaks.Add(p); else pg.NegativeIsoPeaks.Add(p);
            }
            pg.UpdateMonoMassAndIsotopeIntensities(Polarity.Positive);
            foreach (var l in File.ReadAllLines(Path.Combine(DiffDir, "monomass_cpp_result.txt")))
            {
                if (string.IsNullOrWhiteSpace(l)) continue;
                var t = l.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (t[0] == "mono") AssertClose(pg.MonoisotopicMass, t[1], "monoMass");
                else if (t[0] == "intensity") AssertClose(pg.Intensity, t[1], "intensity");
                else AssertClose(pg.PerIsotopeInt[int.Parse(t[1], CultureInfo.InvariantCulture)], t[2], $"perIso[{t[1]}]");
            }
        }

        [Test]
        public void TrimAndNormalize_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "trim_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run trim_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            double[] raw = ParseRow(lines[1]);

            double[] b = MetaFlashDeconAveragine.TrimAndNormalize(raw, out int apex, out int left, out int right);

            var cpp = File.ReadAllLines(Path.Combine(DiffDir, "trim_cpp_result.txt"));
            var hdr = cpp[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            Assert.That(apex, Is.EqualTo(int.Parse(hdr[0], CultureInfo.InvariantCulture)), "apexIndex");
            Assert.That(left, Is.EqualTo(int.Parse(hdr[1], CultureInfo.InvariantCulture)), "leftCountFromApex");
            Assert.That(right, Is.EqualTo(int.Parse(hdr[2], CultureInfo.InvariantCulture)), "rightCountFromApex");
            Assert.That(b.Length, Is.EqualTo(int.Parse(hdr[3], CultureInfo.InvariantCulture)), "b length");
            double[] cppB = ParseRow(cpp[1]);
            for (int i = 0; i < b.Length; i++)
                Assert.That(b[i], Is.EqualTo(cppB[i]).Within(1e-7), $"b[{i}]");
        }

        [Test]
        public void PerChargeCos_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "percc_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run percc_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int minC = int.Parse(h[0], CultureInfo.InvariantCulture);
            int maxC = int.Parse(h[1], CultureInfo.InvariantCulture);
            int curSize = int.Parse(h[2], CultureInfo.InvariantCulture);
            double[] b = ParseRow(lines[1]);

            var pg = new MetaFlashDeconPeakGroup
            {
                MinAbsCharge = minC, MaxAbsCharge = maxC, MinNegativeIsotopeIndex = -1,
                // PerIsotopeInt.Length + MinNegativeIsotopeIndex must equal curSize
                PerIsotopeInt = new double[curSize - (-1)],
            };
            for (int i = 2; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                int z = int.Parse(t[0], CultureInfo.InvariantCulture);
                int iso = int.Parse(t[1], CultureInfo.InvariantCulture);
                double inten = double.Parse(t[2], CultureInfo.InvariantCulture);
                pg.SignalPeaks.Add(new MetaFlashDeconAlgorithm.LogMzPeak(0.0, inten, 0.0, z, iso));
            }

            pg.UpdatePerChargeCos(b);

            foreach (var l in File.ReadAllLines(Path.Combine(DiffDir, "percc_cpp_result.txt")))
            {
                if (string.IsNullOrWhiteSpace(l)) continue;
                var t = l.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                int z = int.Parse(t[0], CultureInfo.InvariantCulture);
                double cpp = double.Parse(t[1], CultureInfo.InvariantCulture);
                Assert.That(pg.PerChargeCos[z], Is.EqualTo(cpp).Within(1e-4), $"perChargeCos z{z}");
            }
        }

        [Test]
        public void IsotopeCosineOffset_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "isocos_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run isocos_cpp.exe first)");
            var lines = File.ReadAllLines(inputPath);
            var h = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int apexIndex = int.Parse(h[0], CultureInfo.InvariantCulture);
            int isoIntShift = int.Parse(h[1], CultureInfo.InvariantCulture);
            int windowWidth = int.Parse(h[2], CultureInfo.InvariantCulture);
            int minIsoSize = int.Parse(h[3], CultureInfo.InvariantCulture);
            double[] per = ParseRow(lines[1]);
            double[] b = ParseRow(lines[2]);

            double cos = MetaFlashDeconPeakGroup.GetIsotopeCosineAndDetermineIsotopeIndex(
                per, b, apexIndex, isoIntShift, windowWidth, minIsoSize, out int offset);

            var cpp = File.ReadAllText(Path.Combine(DiffDir, "isocos_cpp_result.txt")).Split(' ', StringSplitOptions.RemoveEmptyEntries);
            Assert.That(cos, Is.EqualTo(double.Parse(cpp[0], CultureInfo.InvariantCulture)).Within(1e-4), "isoCos");
            Assert.That(offset, Is.EqualTo(int.Parse(cpp[1], CultureInfo.InvariantCulture)), "offset");
        }

        [Test]
        public void Qscore_MatchesOpenMS()
        {
            string inputPath = Path.Combine(DiffDir, "qscore_input.txt");
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run qscore_cpp.exe first)");
            var inLines = File.ReadAllLines(inputPath);
            var cppLines = File.ReadAllLines(Path.Combine(DiffDir, "qscore_cpp_result.txt"));
            int n = 0;
            for (int i = 0; i < inLines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(inLines[i])) continue;
                var t = inLines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                double cos = double.Parse(t[0], CultureInfo.InvariantCulture);
                double snr = double.Parse(t[1], CultureInfo.InvariantCulture);
                double cs = MetaFlashDeconPeakGroup.ComputeQscore(cos, snr);
                double cpp = double.Parse(cppLines[n].Trim(), CultureInfo.InvariantCulture);
                Assert.That(cs, Is.EqualTo(cpp).Within(1e-5), $"qscore(cos={cos},snr={snr})");
                n++;
            }
        }

        private static void AssertClose(double cs, string cppStr, string what)
        {
            double cpp = double.Parse(cppStr, CultureInfo.InvariantCulture);
            // OpenMS computes noise/SNR in float; C# in double -> compare with a magnitude-relative tol.
            double tol = Math.Abs(cpp) * 1e-4 + 1e-2;
            Assert.That(cs, Is.EqualTo(cpp).Within(tol), $"{what}: C# {cs} != OpenMS {cpp}");
        }

        private static double[] ParseRow(string line)
        {
            var t = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            var v = new double[t.Length];
            for (int i = 0; i < t.Length; i++) v[i] = double.Parse(t[i], CultureInfo.InvariantCulture);
            return v;
        }

        private static void RunCase(string inputName, string cppResultName, string csResultName)
        {
            string inputPath = Path.Combine(DiffDir, inputName);
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run noise_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var header = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int absCharge = int.Parse(header[0], CultureInfo.InvariantCulture);
            double isoDa = double.Parse(header[1], CultureInfo.InvariantCulture);

            var signal = new List<(double mz, double intensity)>();
            var noisy = new List<(double mz, double intensity)>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                double mz = double.Parse(t[1], CultureInfo.InvariantCulture);
                double inten = double.Parse(t[2], CultureInfo.InvariantCulture);
                if (t[0] == "S") signal.Add((mz, inten));
                else noisy.Add((mz, inten));
            }

            double cs = MetaFlashDeconAlgorithm.ComputeNoisePeakPower(noisy, signal, absCharge, isoDa);
            TestContext.Progress.WriteLine($"CS  {inputName} getNoisePeakPower = {cs:F8}");
            File.WriteAllText(Path.Combine(DiffDir, csResultName),
                $"{cs:F8}" + Environment.NewLine);

            // Compare against the C++ golden value if present.
            string cppResultPath = Path.Combine(DiffDir, cppResultName);
            if (File.Exists(cppResultPath))
            {
                double cpp = double.Parse(File.ReadAllText(cppResultPath).Trim(), CultureInfo.InvariantCulture);
                Assert.That(cs, Is.EqualTo(cpp).Within(1e-3),
                    $"C# ({cs}) != OpenMS ({cpp}) for getNoisePeakPower");
            }
        }
    }
}
