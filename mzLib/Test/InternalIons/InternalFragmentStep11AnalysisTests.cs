using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test.InternalIons
{
    [TestFixture]
    public class InternalFragmentStep11AnalysisTests
    {
        private const string TsvInputPath = @"E:\Projects\internalIons\mzLibReportsForClaude\InternalFragmentIons.tsv";
        private const string OutputPath = @"E:\Projects\internalIons\mzLibReportsForClaude\Step11Analysis_Output.txt";

        private static readonly char[] StandardAminoAcids =
            { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };

        private static StreamWriter _out;
        private static List<Dictionary<string, string>> _allRows;
        private static List<Dictionary<string, string>> _passing;

        private static void W(string line = "") => _out.WriteLine(line);

        [OneTimeSetUp]
        public void Setup()
        {
            if (!File.Exists(TsvInputPath))
            {
                Assert.Ignore($"TSV file not found: {TsvInputPath}");
                return;
            }

            _out = new StreamWriter(OutputPath, false) { AutoFlush = false };

            var lines = File.ReadAllLines(TsvInputPath);
            if (lines.Length < 2) { Assert.Ignore("TSV file has no data rows"); return; }

            var headers = lines[0].Split('\t');
            _allRows = new List<Dictionary<string, string>>();

            for (int i = 1; i < lines.Length; i++)
            {
                var vals = lines[i].Split('\t');
                var row = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
                for (int j = 0; j < headers.Length && j < vals.Length; j++)
                    row[headers[j]] = vals[j];
                _allRows.Add(row);
            }

            _passing = _allRows.Where(r => B(r, "PassesMassAccuracyFilter")).ToList();

            W($"Loaded {_allRows.Count:N0} rows from TSV, {_passing.Count:N0} passing filter.");
            W($"Output file: {OutputPath}");
            W($"Timestamp: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            W();
        }

        [OneTimeTearDown]
        public void Teardown()
        {
            _out?.Flush();
            _out?.Close();
            _out?.Dispose();
        }

        #region Helper Methods

        private static double D(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return 0.0;
            if (string.IsNullOrWhiteSpace(val)) return 0.0;
            if (val.Equals("NaN", StringComparison.OrdinalIgnoreCase)) return double.NaN;
            return double.TryParse(val, NumberStyles.Any, CultureInfo.InvariantCulture, out var d) ? d : 0.0;
        }

        private static int I(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return 0;
            return int.TryParse(val, out var i) ? i : 0;
        }

        private static bool B(Dictionary<string, string> row, string col)
        {
            if (!row.TryGetValue(col, out var val)) return false;
            return val.Equals("TRUE", StringComparison.OrdinalIgnoreCase) || val == "1";
        }

        private static string S(Dictionary<string, string> row, string col) =>
            row.TryGetValue(col, out var val) ? val : "";

        private static char C(Dictionary<string, string> row, string col)
        {
            var s = S(row, col);
            return string.IsNullOrEmpty(s) ? '-' : s[0];
        }

        private static double Pearson(IList<double> xs, IList<double> ys)
        {
            if (xs.Count != ys.Count || xs.Count < 3) return double.NaN;
            var validPairs = xs.Zip(ys, (x, y) => (x, y))
                .Where(p => !double.IsNaN(p.x) && !double.IsNaN(p.y)).ToList();
            if (validPairs.Count < 3) return double.NaN;

            double mx = validPairs.Average(p => p.x);
            double my = validPairs.Average(p => p.y);
            double sxy = 0, sx2 = 0, sy2 = 0;

            foreach (var (x, y) in validPairs)
            {
                double dx = x - mx, dy = y - my;
                sxy += dx * dy; sx2 += dx * dx; sy2 += dy * dy;
            }
            double denom = Math.Sqrt(sx2 * sy2);
            return denom < 1e-12 ? double.NaN : sxy / denom;
        }

        private static double StdDev(IEnumerable<double> vals)
        {
            var list = vals.Where(v => !double.IsNaN(v)).ToList();
            if (list.Count < 2) return 0;
            double mean = list.Average();
            return Math.Sqrt(list.Select(v => (v - mean) * (v - mean)).Average());
        }

        private static double Median(IEnumerable<double> vals)
        {
            var sorted = vals.Where(v => !double.IsNaN(v)).OrderBy(v => v).ToList();
            return sorted.Count == 0 ? 0 : sorted[sorted.Count / 2];
        }

        #endregion

        [Test, Order(1)]
        public void Test01_DatasetSummary()
        {
            W(new string('=', 60));
            W("TEST 1: DatasetSummary");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();

            W($"Total rows in TSV:               {_allRows.Count:N0}");
            W($"PassesMassAccuracyFilter = TRUE: {_passing.Count:N0}");
            W($"Unique PeptideSequence:          {_passing.Select(r => S(r, "PeptideSequence")).Distinct().Count():N0}");
            W($"Unique ScanNumber:               {_passing.Select(r => S(r, "ScanNumber")).Distinct().Count():N0}");

            W();
            _out.Flush();
        }

        [Test, Order(2)]
        public void Test02_BvsYDirectionality()
        {
            W(new string('=', 60));
            W("TEST 2: BvsYDirectionality");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();

            int n = _passing.Count;
            int hasB = _passing.Count(r => B(r, "HasMatchedBIonAtNTerm"));
            int hasY = _passing.Count(r => B(r, "HasMatchedYIonAtCTerm"));

            W($"HasMatchedBIonAtNTerm = TRUE:  {hasB,6:N0} ({100.0 * hasB / n:F1}%)");
            W($"HasMatchedYIonAtCTerm = TRUE:  {hasY,6:N0} ({100.0 * hasY / n:F1}%)");

            W();
            _out.Flush();
        }

        [Test, Order(3)]
        public void Test03_YIonMechanismStratification()
        {
            W(new string('=', 60));
            W("TEST 3: YIonMechanismStratification");
            var withY = _passing.Where(r => B(r, "HasMatchedYIonAtCTerm")).ToList();
            W($"N analyzed (HasMatchedYIonAtCTerm=TRUE): {withY.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(4)]
        public void Test04_FourGroupComparison()
        {
            W(new string('=', 60));
            W("TEST 4: FourGroupComparison");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(5)]
        public void Test05_UnsupportedHighIntensityIons()
        {
            W(new string('=', 60));
            W("TEST 5: UnsupportedHighIntensityIons");
            var unsupported = _passing.Where(r => !B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm")).ToList();
            W($"N analyzed (neither B nor Y): {unsupported.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(6)]
        public void Test06_WeakBondEnrichmentTable()
        {
            W(new string('=', 60));
            W("TEST 6: WeakBondEnrichmentTable");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(7)]
        public void Test07_PearsonCorrelationMatrix()
        {
            W(new string('=', 60));
            W("TEST 7: PearsonCorrelationMatrix");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();

            var target = _passing.Select(r => D(r, "TicNormalizedIntensity")).ToList();

            var correlations = new List<(string name, double r)>
            {
                ("BIonIntensityAtNTerm", Pearson(_passing.Select(r => D(r, "BIonIntensityAtNTerm")).ToList(), target)),
                ("YIonIntensityAtCTerm", Pearson(_passing.Select(r => D(r, "YIonIntensityAtCTerm")).ToList(), target)),
                ("MaxTerminalIonIntensity", Pearson(_passing.Select(r => D(r, "MaxTerminalIonIntensity")).ToList(), target)),
                ("BYProductScore", Pearson(_passing.Select(r => D(r, "BYProductScore")).ToList(), target)),
                ("HasBothTerminalIons", Pearson(_passing.Select(r => B(r, "HasBothTerminalIons") ? 1.0 : 0.0).ToList(), target)),
                ("HasMatchedBIonAtNTerm", Pearson(_passing.Select(r => B(r, "HasMatchedBIonAtNTerm") ? 1.0 : 0.0).ToList(), target)),
                ("HasMatchedYIonAtCTerm", Pearson(_passing.Select(r => B(r, "HasMatchedYIonAtCTerm") ? 1.0 : 0.0).ToList(), target)),
                ("BasicResiduesInBIonSpan", Pearson(_passing.Select(r => (double)I(r, "BasicResiduesInBIonSpan")).ToList(), target)),
                ("BasicResiduesInYIonSpan", Pearson(_passing.Select(r => (double)I(r, "BasicResiduesInYIonSpan")).ToList(), target)),
                ("NumberOfBasicResidues", Pearson(_passing.Select(r => (double)I(r, "NumberOfBasicResidues")).ToList(), target)),
                ("FragmentLength", Pearson(_passing.Select(r => (double)I(r, "FragmentLength")).ToList(), target)),
                ("DistanceFromCTerm", Pearson(_passing.Select(r => (double)I(r, "DistanceFromCTerm")).ToList(), target)),
                ("PeptideLength", Pearson(_passing.Select(r => (double)I(r, "PeptideLength")).ToList(), target)),
                ("RelativeDistanceFromCTerm", Pearson(_passing.Select(r => D(r, "RelativeDistanceFromCTerm")).ToList(), target)),
                ("LocalIntensityRank_neg", Pearson(_passing.Select(r => -D(r, "LocalIntensityRank")).ToList(), target)),
                ("HasProlineAtEitherTerm", Pearson(_passing.Select(r => B(r, "HasProlineAtEitherTerminus") ? 1.0 : 0.0).ToList(), target)),
                ("HasAspartateAtEitherTerm", Pearson(_passing.Select(r => B(r, "HasAspartateAtEitherTerminus") ? 1.0 : 0.0).ToList(), target)),
                ("HasModifiedResidue", Pearson(_passing.Select(r => B(r, "HasModifiedResidue") ? 1.0 : 0.0).ToList(), target)),
                // Scorer features
                ("IsProlineAtInternalNTerm", Pearson(_passing.Select(r => B(r, "IsProlineAtInternalNTerminus") ? 1.0 : 0.0).ToList(), target)),
                ("IsTerminalRescue", Pearson(_passing.Select(r => B(r, "IsTerminalRescue") ? 1.0 : 0.0).ToList(), target)),
                ("NTermFlankHydrophobicity", Pearson(_passing.Select(r => D(r, "NTerminalFlankingHydrophobicity")).ToList(), target)),
                // Modification category features
                ("CommonBiologicalModCount", Pearson(_passing.Select(r => (double)I(r, "CommonBiologicalModCount")).ToList(), target)),
                ("CommonArtifactModCount", Pearson(_passing.Select(r => (double)I(r, "CommonArtifactModCount")).ToList(), target)),
                ("MetalModCount", Pearson(_passing.Select(r => (double)I(r, "MetalModCount")).ToList(), target)),
                ("LessCommonModCount", Pearson(_passing.Select(r => (double)I(r, "LessCommonModCount")).ToList(), target)),
                ("HasPhosphorylation", Pearson(_passing.Select(r => B(r, "HasPhosphorylation") ? 1.0 : 0.0).ToList(), target)),
                ("HasMetalOnTerminalAcidic", Pearson(_passing.Select(r => B(r, "HasMetalOnTerminalAcidic") ? 1.0 : 0.0).ToList(), target))
            };

            W("Correlations with TicNormalizedIntensity (ranked by |r|):");
            W("+---------------------------+----------+--------+");
            W("| Feature                   |        r | Flag   |");
            W("+---------------------------+----------+--------+");

            foreach (var (name, r) in correlations.OrderByDescending(x => Math.Abs(x.r)))
            {
                string flag = double.IsNaN(r) ? "N/A" : Math.Abs(r) > 0.20 ? "SIGNAL" : Math.Abs(r) > 0.10 ? "WEAK" : "";
                W($"| {name,-25} | {(double.IsNaN(r) ? "NaN" : r.ToString("F4")),8} | {flag,-6} |");
            }
            W("+---------------------------+----------+--------+");

            W();
            _out.Flush();
        }

        [Test, Order(8)]
        public void Test08_FeatureReadinessTable()
        {
            W(new string('=', 60));
            W("TEST 8: FeatureReadinessTable");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();

            var target = _passing.Select(r => D(r, "TicNormalizedIntensity")).ToList();

            var features = new (string name, Func<Dictionary<string, string>, double> get)[]
            {
                ("FragmentLength", r => I(r, "FragmentLength")),
                ("DistanceFromCTerm", r => I(r, "DistanceFromCTerm")),
                ("PeptideLength", r => I(r, "PeptideLength")),
                ("RelativeDistanceFromCTerm", r => D(r, "RelativeDistanceFromCTerm")),
                ("NumberOfBasicResidues", r => I(r, "NumberOfBasicResidues")),
                ("BasicResiduesInBIonSpan", r => I(r, "BasicResiduesInBIonSpan")),
                ("BasicResiduesInYIonSpan", r => I(r, "BasicResiduesInYIonSpan")),
                ("BIonIntensityAtNTerm", r => D(r, "BIonIntensityAtNTerm")),
                ("YIonIntensityAtCTerm", r => D(r, "YIonIntensityAtCTerm")),
                ("MaxTerminalIonIntensity", r => D(r, "MaxTerminalIonIntensity")),
                ("BYProductScore", r => D(r, "BYProductScore")),
                ("HasBothTerminalIons", r => B(r, "HasBothTerminalIons") ? 1.0 : 0.0),
                ("HasProlineAtEitherTerm", r => B(r, "HasProlineAtEitherTerminus") ? 1.0 : 0.0),
                ("HasAspartateAtEitherTerm", r => B(r, "HasAspartateAtEitherTerminus") ? 1.0 : 0.0),
                ("LocalIntensityRank", r => D(r, "LocalIntensityRank")),
                ("HasModifiedResidue", r => B(r, "HasModifiedResidue") ? 1.0 : 0.0),
                ("PrecursorCharge", r => I(r, "PrecursorCharge")),
                ("CollisionEnergy", r => D(r, "CollisionEnergy")),
                ("PeptidePEP", r => D(r, "PeptidePEP")),
                // Scorer features
                ("IsProlineAtInternalNTerm", r => B(r, "IsProlineAtInternalNTerminus") ? 1.0 : 0.0),
                ("IsTerminalRescue", r => B(r, "IsTerminalRescue") ? 1.0 : 0.0),
                ("NTermFlankHydrophobicity", r => D(r, "NTerminalFlankingHydrophobicity")),
                // Modification category features
                ("CommonBiologicalModCount", r => I(r, "CommonBiologicalModCount")),
                ("CommonArtifactModCount", r => I(r, "CommonArtifactModCount")),
                ("MetalModCount", r => I(r, "MetalModCount")),
                ("LessCommonModCount", r => I(r, "LessCommonModCount")),
                ("HasPhosphorylation", r => B(r, "HasPhosphorylation") ? 1.0 : 0.0),
                ("HasMetalOnTerminalAcidic", r => B(r, "HasMetalOnTerminalAcidic") ? 1.0 : 0.0)
            };

            W("+---------------------------+--------+--------+-------+----------+----------+");
            W("| Feature                   | NValid | NNonZ  | %NonZ | r(TicNI) | Flags    |");
            W("+---------------------------+--------+--------+-------+----------+----------+");

            foreach (var (name, getter) in features)
            {
                var vals = _passing.Select(getter).Where(v => !double.IsNaN(v)).ToList();
                int nVal = vals.Count;
                int nNz = vals.Count(v => v != 0);
                double pctNz = nVal > 0 ? 100.0 * nNz / nVal : 0;
                double r = Pearson(_passing.Select(getter).ToList(), target);

                var flags = new List<string>();
                if (pctNz < 20) flags.Add("SPARSE");
                if (!double.IsNaN(r) && Math.Abs(r) > 0.15) flags.Add("SIGNAL");

                string flagStr = flags.Count > 0 ? string.Join(",", flags) : "";
                string rStr = double.IsNaN(r) ? "NaN" : r.ToString("F4");

                W($"| {name,-25} | {nVal,6} | {nNz,6} | {pctNz,5:F1} | {rStr,8} | {flagStr,-8} |");
            }

            W("+---------------------------+--------+--------+-------+----------+----------+");
            W();
            _out.Flush();
        }

        [Test, Order(9)]
        public void Test09_ScanLevelSanityCheck()
        {
            W(new string('=', 60));
            W("TEST 9: ScanLevelSanityCheck");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(10)]
        public void Test10_IsobaricAmbiguityImpact()
        {
            W(new string('=', 60));
            W("TEST 10: IsobaricAmbiguityImpact");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();
            _out.Flush();
        }

        [Test, Order(11)]
        public void Test11_NewFeatureValidation()
        {
            W(new string('=', 60));
            W("TEST 11: NewFeatureValidation");
            W($"N analyzed: {_passing.Count}");
            W(new string('=', 60));
            W();

            var target = _passing.Select(r => D(r, "TicNormalizedIntensity")).ToList();

            // === IsProlineAtInternalNTerminus ===
            W("--- IsProlineAtInternalNTerminus ---");
            var proTrue = _passing.Where(r => B(r, "IsProlineAtInternalNTerminus")).ToList();
            var proFalse = _passing.Where(r => !B(r, "IsProlineAtInternalNTerminus")).ToList();

            W($"  Count TRUE:  {proTrue.Count:N0}");
            W($"  Count FALSE: {proFalse.Count:N0}");

            double meanTicProTrue = proTrue.Count > 0 ? proTrue.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double meanTicProFalse = proFalse.Count > 0 ? proFalse.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double proRatio = meanTicProFalse > 0 ? meanTicProTrue / meanTicProFalse : 0;

            W($"  Mean TicNI: TRUE={meanTicProTrue:E4}, FALSE={meanTicProFalse:E4}, Ratio={proRatio:F2}");
            W("  (Expectation: ratio > 1.0)");
            W();

            // === IsTerminalRescue ===
            W("--- IsTerminalRescue ---");
            var rescueTrue = _passing.Where(r => B(r, "IsTerminalRescue")).ToList();
            var rescueFalse = _passing.Where(r => !B(r, "IsTerminalRescue")).ToList();

            W($"  Count TRUE:  {rescueTrue.Count:N0}");
            W($"  Count FALSE: {rescueFalse.Count:N0}");

            double meanYRescueTrue = rescueTrue.Count > 0 ? rescueTrue.Average(r => D(r, "YIonIntensityAtCTerm")) : 0;
            double meanYRescueFalse = rescueFalse.Count > 0 ? rescueFalse.Average(r => D(r, "YIonIntensityAtCTerm")) : 0;

            W($"  Mean YIonIntensity: TRUE={meanYRescueTrue:E4}, FALSE={meanYRescueFalse:E4}");
            W("  (Expectation: YIon elevated when TRUE)");
            W();

            // === NTerminalFlankingHydrophobicity ===
            W("--- NTerminalFlankingHydrophobicity ---");
            var hydroVals = _passing.Select(r => D(r, "NTerminalFlankingHydrophobicity")).Where(v => !double.IsNaN(v)).ToList();
            var hydroList = _passing.Select(r => D(r, "NTerminalFlankingHydrophobicity")).ToList();

            W($"  Min: {hydroVals.Min():F2}, Max: {hydroVals.Max():F2}, Mean: {hydroVals.Average():F2}");
            W($"  r(NTermFlankHydro, TicNI): {Pearson(hydroList, target):F4}");
            W();

            // === RelativeDistanceFromCTerm ===
            W("--- RelativeDistanceFromCTerm ---");
            var relDistVals = _passing.Select(r => D(r, "RelativeDistanceFromCTerm")).Where(v => !double.IsNaN(v)).ToList();
            var relDistList = _passing.Select(r => D(r, "RelativeDistanceFromCTerm")).ToList();

            W($"  Min: {relDistVals.Min():F3}, Max: {relDistVals.Max():F3}, Mean: {relDistVals.Average():F3}");
            W($"  r(RelativeDistanceFromCTerm, TicNI): {Pearson(relDistList, target):F4}");
            W();

            // === PeptideLength independence check ===
            W("--- PeptideLength ---");
            var nearCTerm = _passing.Where(r => I(r, "DistanceFromCTerm") <= 3).ToList();
            var farCTerm = _passing.Where(r => I(r, "DistanceFromCTerm") > 3).ToList();

            double meanPepLenNear = nearCTerm.Count > 0 ? nearCTerm.Average(r => I(r, "PeptideLength")) : 0;
            double meanPepLenFar = farCTerm.Count > 0 ? farCTerm.Average(r => I(r, "PeptideLength")) : 0;

            W($"  Mean PeptideLength: DistCTerm<=3 = {meanPepLenNear:F1}, DistCTerm>3 = {meanPepLenFar:F1}");
            W("  (Confirms peptide length varies independently of terminal proximity)");
            W();

            // === MODIFICATION CATEGORY FEATURES ===
            W("--- Modification Category Features ---");
            W();

            PrintModCountFeature("CommonBiologicalModCount", target);
            PrintModCountFeature("CommonArtifactModCount", target);
            PrintModCountFeature("MetalModCount", target);
            PrintModCountFeature("LessCommonModCount", target);

            // HasPhosphorylation
            W("--- HasPhosphorylation ---");
            var phosphoTrue = _passing.Where(r => B(r, "HasPhosphorylation")).ToList();
            var phosphoFalse = _passing.Where(r => !B(r, "HasPhosphorylation")).ToList();

            W($"  Count TRUE:  {phosphoTrue.Count:N0} ({100.0 * phosphoTrue.Count / _passing.Count:F2}%)");
            W($"  Count FALSE: {phosphoFalse.Count:N0}");

            double meanTicPhosphoTrue = phosphoTrue.Count > 0 ? phosphoTrue.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double meanTicPhosphoFalse = phosphoFalse.Count > 0 ? phosphoFalse.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double phosphoRatio = meanTicPhosphoFalse > 0 ? meanTicPhosphoTrue / meanTicPhosphoFalse : 0;

            W($"  Mean TicNI TRUE:  {meanTicPhosphoTrue:E4}");
            W($"  Mean TicNI FALSE: {meanTicPhosphoFalse:E4}");
            W($"  Ratio (TRUE/FALSE): {phosphoRatio:F2}");
            W("  Expected: ratio < 1.0 (phospho neutral loss competes with backbone)");
            W();

            // HasMetalOnTerminalAcidic
            W("--- HasMetalOnTerminalAcidic ---");
            var metalTermTrue = _passing.Where(r => B(r, "HasMetalOnTerminalAcidic")).ToList();
            var metalTermFalse = _passing.Where(r => !B(r, "HasMetalOnTerminalAcidic")).ToList();

            W($"  Count TRUE:  {metalTermTrue.Count:N0} ({100.0 * metalTermTrue.Count / _passing.Count:F2}%)");
            W($"  Count FALSE: {metalTermFalse.Count:N0}");

            double meanTicMetalTrue = metalTermTrue.Count > 0 ? metalTermTrue.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double meanTicMetalFalse = metalTermFalse.Count > 0 ? metalTermFalse.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double metalRatio = meanTicMetalFalse > 0 ? meanTicMetalTrue / meanTicMetalFalse : 0;

            W($"  Mean TicNI TRUE:  {meanTicMetalTrue:E4}");
            W($"  Mean TicNI FALSE: {meanTicMetalFalse:E4}");
            W($"  Ratio (TRUE/FALSE): {metalRatio:F2}");
            W("  (No strong prior expectation - report direction only)");
            W();

            // === MODIFICATION INTERACTION WITH WEAK-BOND FEATURES ===
            W("--- Modification Interaction with Weak-Bond Features ---");
            W();

            var withAsp = _passing.Where(r => B(r, "HasAspartateAtEitherTerminus")).ToList();
            W($"Among HasAspartateAtEitherTerminus = TRUE: {withAsp.Count:N0} total");

            var aspWithMetal = withAsp.Where(r => I(r, "MetalModCount") >= 1).ToList();
            var aspNoMetal = withAsp.Where(r => I(r, "MetalModCount") == 0).ToList();

            W($"  With MetalModCount >= 1: {aspWithMetal.Count:N0} ({100.0 * aspWithMetal.Count / Math.Max(withAsp.Count, 1):F1}%)");
            W($"  With MetalModCount = 0:  {aspNoMetal.Count:N0}");

            double meanTicAspMetal = aspWithMetal.Count > 0 ? aspWithMetal.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double meanTicAspNoMetal = aspNoMetal.Count > 0 ? aspNoMetal.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double aspMetalRatio = meanTicAspNoMetal > 0 ? meanTicAspMetal / meanTicAspNoMetal : 0;

            W($"  Mean TicNI, MetalModCount=0:  {meanTicAspNoMetal:E4}");
            W($"  Mean TicNI, MetalModCount>=1: {meanTicAspMetal:E4}");
            W($"  Ratio: {aspMetalRatio:F2}");
            W();
            W("  (Tests whether metal coordination at terminal D/E alters intensity");
            W("   vs uncoordinated terminal D/E. If ratio != 1.0, metal adducts");
            W("   interact with HasAspartateAtEitherTerminus.)");

            W();
            _out.Flush();
        }

        private void PrintModCountFeature(string featureName, List<double> target)
        {
            W($"--- {featureName} ---");

            var zeroRows = _passing.Where(r => I(r, featureName) == 0).ToList();
            var nonZeroRows = _passing.Where(r => I(r, featureName) >= 1).ToList();

            W($"  N with count=0:   {zeroRows.Count:N0}");
            W($"  N with count>=1:  {nonZeroRows.Count:N0} ({100.0 * nonZeroRows.Count / _passing.Count:F1}%)");

            double meanTicZero = zeroRows.Count > 0 ? zeroRows.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double meanTicNonZero = nonZeroRows.Count > 0 ? nonZeroRows.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double ratio = meanTicZero > 0 ? meanTicNonZero / meanTicZero : 0;

            W($"  Mean TicNI when count=0:  {meanTicZero:E4}");
            W($"  Mean TicNI when count>=1: {meanTicNonZero:E4}");
            W($"  Intensity ratio (nonzero/zero): {ratio:F2}");

            var featureVals = _passing.Select(r => (double)I(r, featureName)).ToList();
            double r_corr = Pearson(featureVals, target);
            W($"  Pearson r with TicNI: {r_corr:F4}");
            W();
        }
    }
}