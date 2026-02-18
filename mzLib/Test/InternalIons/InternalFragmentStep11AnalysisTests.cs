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
            if (lines.Length < 2)
            {
                Assert.Ignore("TSV file has no data rows");
                return;
            }

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

        private static string S(Dictionary<string, string> row, string col)
        {
            return row.TryGetValue(col, out var val) ? val : "";
        }

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
                sxy += dx * dy;
                sx2 += dx * dx;
                sy2 += dy * dy;
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
            if (sorted.Count == 0) return 0;
            return sorted[sorted.Count / 2];
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

            var ticNI = _passing.Select(r => D(r, "TicNormalizedIntensity")).Where(v => !double.IsNaN(v)).ToList();
            W("TicNormalizedIntensity (passing):");
            W($"  Min:    {ticNI.Min():E4}");
            W($"  Max:    {ticNI.Max():E4}");
            W($"  Mean:   {ticNI.Average():E4}");
            W($"  StdDev: {StdDev(ticNI):E4}");
            W();

            var tics = _passing.Select(r => D(r, "TotalIonCurrent")).Where(v => v > 0).Distinct().ToList();
            if (tics.Count > 0)
            {
                W($"TotalIonCurrent: Min={tics.Min():E4}, Max={tics.Max():E4}");
                W();
            }

            W("FragmentLength distribution:");
            for (int len = 3; len <= 24; len++)
            {
                int count = _passing.Count(r => I(r, "FragmentLength") == len);
                if (count > 0) W($"  Len {len,2}: {count,6:N0} fragments");
            }
            int over25 = _passing.Count(r => I(r, "FragmentLength") >= 25);
            if (over25 > 0) W($"  Len 25+: {over25,6:N0} fragments");
            W();

            W($"IsIsobaricAmbiguous = TRUE:  {_passing.Count(r => B(r, "IsIsobaricAmbiguous")):N0}");
            W($"HasModifiedResidue = TRUE:   {_passing.Count(r => B(r, "HasModifiedResidue")):N0}");

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
            int hasBoth = _passing.Count(r => B(r, "HasBothTerminalIons"));
            int hasNeither = _passing.Count(r => !B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm"));

            W($"HasMatchedBIonAtNTerm = TRUE:  {hasB,6:N0} ({100.0 * hasB / n:F1}%)");
            W($"HasMatchedYIonAtCTerm = TRUE:  {hasY,6:N0} ({100.0 * hasY / n:F1}%)");
            W($"HasBothTerminalIons = TRUE:    {hasBoth,6:N0} ({100.0 * hasBoth / n:F1}%)");
            W($"Neither matched:               {hasNeither,6:N0} ({100.0 * hasNeither / n:F1}%)");
            W();

            var bVals = _passing.Select(r => D(r, "BIonIntensityAtNTerm")).Where(v => v > 0).ToList();
            var yVals = _passing.Select(r => D(r, "YIonIntensityAtCTerm")).Where(v => v > 0).ToList();

            double meanB = bVals.Count > 0 ? bVals.Average() : 0;
            double meanY = yVals.Count > 0 ? yVals.Average() : 0;

            W($"BIonIntensityAtNTerm (B>0, n={bVals.Count:N0}):");
            W($"  Mean: {meanB:E4}, StdDev: {StdDev(bVals):E4}, Max: {(bVals.Count > 0 ? bVals.Max() : 0):E4}");
            W($"YIonIntensityAtCTerm (Y>0, n={yVals.Count:N0}):");
            W($"  Mean: {meanY:E4}, StdDev: {StdDev(yVals):E4}, Max: {(yVals.Count > 0 ? yVals.Max() : 0):E4}");
            W();

            if (meanB > meanY && meanB > 0)
                W($"CONFIRMED: B-ions dominate (meanB={meanB:E4} > meanY={meanY:E4})");
            else
                W($"UNEXPECTED: Y-ions match or exceed B-ions (meanB={meanB:E4}, meanY={meanY:E4})");

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

            var s1 = withY.Where(r => I(r, "DistanceFromCTerm") <= 3).ToList();
            var s2 = withY.Where(r => I(r, "DistanceFromCTerm") > 3 && I(r, "BasicResiduesInYIonSpan") >= 2).ToList();
            var s3 = withY.Where(r => I(r, "DistanceFromCTerm") > 3 && I(r, "BasicResiduesInYIonSpan") <= 1).ToList();

            void PrintStratum(string name, string desc, List<Dictionary<string, string>> list)
            {
                W($"{name}: {list.Count:N0} ions - {desc}");
                if (list.Count == 0) { W(); return; }

                var yVals = list.Select(r => D(r, "YIonIntensityAtCTerm")).ToList();
                var ticVals = list.Select(r => D(r, "TicNormalizedIntensity")).ToList();

                W($"  MeanYIonTic: {yVals.Average():E4}, StdDev: {StdDev(yVals):E4}");
                W($"  MeanTicNI:   {ticVals.Average():E4}, StdDev: {StdDev(ticVals):E4}");
                W($"  Frac TicNI>0.003: {100.0 * list.Count(r => D(r, "TicNormalizedIntensity") > 0.003) / list.Count:F1}%");
                W($"  Frac TicNI>0.01:  {100.0 * list.Count(r => D(r, "TicNormalizedIntensity") > 0.01) / list.Count:F1}%");
                W();

                W("  Top 5 by TicNI:");
                foreach (var r in list.OrderByDescending(r => D(r, "TicNormalizedIntensity")).Take(5))
                {
                    string pep = S(r, "PeptideSequence");
                    if (pep.Length > 20) pep = pep.Substring(0, 17) + "...";
                    W($"    {S(r, "InternalSequence"),-12} TicNI={D(r, "TicNormalizedIntensity"):E3} YIon={D(r, "YIonIntensityAtCTerm"):E3} Dist={I(r, "DistanceFromCTerm")} BasicY={I(r, "BasicResiduesInYIonSpan")} Pep={pep}");
                }
                W();
            }

            PrintStratum("S1", "C-terminal K/R rescue (DistCTerm <= 3)", s1);
            PrintStratum("S2", "Multiply-basic (DistCTerm > 3, BasicInY >= 2)", s2);
            PrintStratum("S3", "Weak retention (DistCTerm > 3, BasicInY <= 1)", s3);

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

            var g00 = _passing.Where(r => !B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm")).ToList();
            var g10 = _passing.Where(r => B(r, "HasMatchedBIonAtNTerm") && !B(r, "HasMatchedYIonAtCTerm")).ToList();
            var g01 = _passing.Where(r => !B(r, "HasMatchedBIonAtNTerm") && B(r, "HasMatchedYIonAtCTerm")).ToList();
            var g11 = _passing.Where(r => B(r, "HasMatchedBIonAtNTerm") && B(r, "HasMatchedYIonAtCTerm")).ToList();

            W("+--------+--------+-----------+-----------+-----------+-----------+----------+----------+");
            W("| Group  |  Count |   MeanTic |  MedianTic|  StdDevTic|    MaxTic | %>0.003  | %>0.01   |");
            W("+--------+--------+-----------+-----------+-----------+-----------+----------+----------+");

            void PrintGroup(string name, List<Dictionary<string, string>> grp)
            {
                if (grp.Count == 0)
                {
                    W($"| {name,-6} |      0 |       N/A |       N/A |       N/A |       N/A |      N/A |      N/A |");
                    return;
                }
                var vals = grp.Select(r => D(r, "TicNormalizedIntensity")).ToList();
                double mean = vals.Average();
                double median = Median(vals);
                double std = StdDev(vals);
                double max = vals.Max();
                double pct003 = 100.0 * vals.Count(v => v > 0.003) / vals.Count;
                double pct01 = 100.0 * vals.Count(v => v > 0.01) / vals.Count;
                W($"| {name,-6} | {grp.Count,6} | {mean,9:E3} | {median,9:E3} | {std,9:E3} | {max,9:E3} | {pct003,7:F1}% | {pct01,7:F1}% |");
            }

            PrintGroup("G00", g00);
            PrintGroup("G10", g10);
            PrintGroup("G01", g01);
            PrintGroup("G11", g11);
            W("+--------+--------+-----------+-----------+-----------+-----------+----------+----------+");
            W();

            double mean10 = g10.Count > 0 ? g10.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double mean01 = g01.Count > 0 ? g01.Average(r => D(r, "TicNormalizedIntensity")) : 0;
            double mean11 = g11.Count > 0 ? g11.Average(r => D(r, "TicNormalizedIntensity")) : 0;

            W($"B-dominance test: G10({mean10:E3}) vs G01({mean01:E3}) => {(mean10 > mean01 ? "G10 > G01 (B dominates)" : "G01 >= G10")}");
            W($"Both-termini test: G11({mean11:E3}) vs G10({mean10:E3}) => {(mean11 > mean10 ? "G11 > G10 (amplification)" : "G10 >= G11")}");

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

            if (unsupported.Count == 0)
            {
                W("No unsupported ions.");
                W();
                _out.Flush();
                return;
            }

            W("Top 20 by TicNormalizedIntensity:");
            W("+------------------+----------+------+------+------+------+-----+-----+-----+----------------------+");
            W("| Sequence         | TicNI    | NFlnk| CFlnk| NtAA | CtAA | Len | Pro | Asp | Peptide              |");
            W("+------------------+----------+------+------+------+------+-----+-----+-----+----------------------+");

            foreach (var r in unsupported.OrderByDescending(r => D(r, "TicNormalizedIntensity")).Take(20))
            {
                string seq = S(r, "InternalSequence");
                if (seq.Length > 16) seq = seq.Substring(0, 13) + "...";
                string pep = S(r, "PeptideSequence");
                if (pep.Length > 20) pep = pep.Substring(0, 17) + "...";

                W($"| {seq,-16} | {D(r, "TicNormalizedIntensity"),8:E2} | {C(r, "NTerminalFlankingResidue"),4} | {C(r, "CTerminalFlankingResidue"),4} | {C(r, "InternalNTerminalAA"),4} | {C(r, "InternalCTerminalAA"),4} | {I(r, "FragmentLength"),3} | {(B(r, "HasProlineAtEitherTerminus") ? "Y" : "N"),3} | {(B(r, "HasAspartateAtEitherTerminus") ? "Y" : "N"),3} | {pep,-20} |");
            }
            W("+------------------+----------+------+------+------+------+-----+-----+-----+----------------------+");
            W();

            W("Weak-bond hypothesis indicators:");
            W($"  HasAspartateAtEitherTerminus = TRUE: {unsupported.Count(r => B(r, "HasAspartateAtEitherTerminus")):N0} ({100.0 * unsupported.Count(r => B(r, "HasAspartateAtEitherTerminus")) / unsupported.Count:F1}%)");
            W($"  HasProlineAtEitherTerminus = TRUE:   {unsupported.Count(r => B(r, "HasProlineAtEitherTerminus")):N0} ({100.0 * unsupported.Count(r => B(r, "HasProlineAtEitherTerminus")) / unsupported.Count:F1}%)");
            W($"  NTerminalFlankingResidue in D/E/P:   {unsupported.Count(r => "DEP".Contains(C(r, "NTerminalFlankingResidue"))):N0}");
            W($"  CTerminalFlankingResidue in D/E/P:   {unsupported.Count(r => "DEP".Contains(C(r, "CTerminalFlankingResidue"))):N0}");
            W($"  InternalNTerminalAA in D/E/P:        {unsupported.Count(r => "DEP".Contains(C(r, "InternalNTerminalAA"))):N0}");

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

            // Background frequency from all internal sequences
            string allSeq = string.Concat(_passing.Select(r => S(r, "InternalSequence")));
            var bgFreq = StandardAminoAcids.ToDictionary(aa => aa,
                aa => allSeq.Length > 0 ? (double)allSeq.Count(c => c == aa) / allSeq.Length : 0);

            void PrintEnrichment(string posName, Func<Dictionary<string, string>, char> getAA)
            {
                W($"\n{posName}:");
                W("+----+-------+---------+--------+----------+--------+");
                W("| AA | Count | Expected| Enrich | MeanTicNI| Flag   |");
                W("+----+-------+---------+--------+----------+--------+");

                var counts = StandardAminoAcids.ToDictionary(aa => aa, aa => _passing.Count(r => getAA(r) == aa));
                int total = _passing.Count;

                var enrichments = StandardAminoAcids
                    .Select(aa =>
                    {
                        double exp = bgFreq[aa] * total;
                        double enrich = exp > 0 ? counts[aa] / exp : 0;
                        var matching = _passing.Where(r => getAA(r) == aa).ToList();
                        double meanTic = matching.Count > 0 ? matching.Average(r => D(r, "TicNormalizedIntensity")) : 0;
                        string flag = enrich > 1.5 ? "ENRICH" : enrich < 0.67 ? "DEPLET" : "";
                        return (aa, count: counts[aa], exp, enrich, meanTic, flag);
                    })
                    .Where(x => x.count > 0)
                    .OrderByDescending(x => x.enrich)
                    .Take(12);

                foreach (var (aa, count, exp, enrich, meanTic, flag) in enrichments)
                    W($"| {aa}  | {count,5} | {exp,7:F1} | {enrich,6:F2} | {meanTic,8:E2} | {flag,-6} |");

                W("+----+-------+---------+--------+----------+--------+");
            }

            PrintEnrichment("NTerminalFlankingResidue", r => C(r, "NTerminalFlankingResidue"));
            PrintEnrichment("CTerminalFlankingResidue", r => C(r, "CTerminalFlankingResidue"));
            PrintEnrichment("InternalNTerminalAA", r => C(r, "InternalNTerminalAA"));
            PrintEnrichment("InternalCTerminalAA", r => C(r, "InternalCTerminalAA"));

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
                ("LocalIntensityRank_neg", Pearson(_passing.Select(r => -D(r, "LocalIntensityRank")).ToList(), target)),
                ("HasProlineAtEitherTerm", Pearson(_passing.Select(r => B(r, "HasProlineAtEitherTerminus") ? 1.0 : 0.0).ToList(), target)),
                ("HasAspartateAtEitherTerm", Pearson(_passing.Select(r => B(r, "HasAspartateAtEitherTerminus") ? 1.0 : 0.0).ToList(), target)),
                ("HasModifiedResidue", Pearson(_passing.Select(r => B(r, "HasModifiedResidue") ? 1.0 : 0.0).ToList(), target))
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

            // Additional mechanism correlations
            var bSpan = _passing.Select(r => (double)I(r, "BasicResiduesInBIonSpan")).ToList();
            var ySpan = _passing.Select(r => (double)I(r, "BasicResiduesInYIonSpan")).ToList();
            var bInt = _passing.Select(r => D(r, "BIonIntensityAtNTerm")).ToList();
            var yInt = _passing.Select(r => D(r, "YIonIntensityAtCTerm")).ToList();

            W("Mechanism validation correlations:");
            W($"  r(BasicResiduesInBIonSpan, BIonIntensity): {Pearson(bSpan, bInt):F4}");
            W($"  r(BasicResiduesInYIonSpan, YIonIntensity): {Pearson(ySpan, yInt):F4}");
            W($"  r(BIonIntensity, YIonIntensity):           {Pearson(bInt, yInt):F4}");

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
                ("PeptidePEP", r => D(r, "PeptidePEP"))
            };

            W("+---------------------------+--------+--------+-------+----------+----------+----------+----------+----------+--------+");
            W("| Feature                   | NValid | NNonZ  | %NonZ |     Mean |   StdDev |      Min |      Max | r(TicNI) | Flags  |");
            W("+---------------------------+--------+--------+-------+----------+----------+----------+----------+----------+--------+");

            foreach (var (name, getter) in features)
            {
                var vals = _passing.Select(getter).Where(v => !double.IsNaN(v)).ToList();
                int nVal = vals.Count;
                int nNz = vals.Count(v => v != 0);
                double pctNz = nVal > 0 ? 100.0 * nNz / nVal : 0;
                double mean = nVal > 0 ? vals.Average() : 0;
                double std = StdDev(vals);
                double min = nVal > 0 ? vals.Min() : 0;
                double max = nVal > 0 ? vals.Max() : 0;
                double r = Pearson(_passing.Select(getter).ToList(), target);

                var flags = new List<string>();
                if (pctNz < 20) flags.Add("SPARSE");
                if (!double.IsNaN(r) && Math.Abs(r) > 0.15) flags.Add("SIGNAL");

                string flagStr = flags.Count > 0 ? string.Join(",", flags) : "";
                string rStr = double.IsNaN(r) ? "NaN" : r.ToString("F4");

                W($"| {name,-25} | {nVal,6} | {nNz,6} | {pctNz,5:F1} | {mean,8:E2} | {std,8:E2} | {min,8:E2} | {max,8:E2} | {rStr,8} | {flagStr,-6} |");
            }

            // Target row
            var targetVals = target.Where(v => !double.IsNaN(v)).ToList();
            W($"| {"TicNormalizedIntensity",-25} | {targetVals.Count,6} | {targetVals.Count(v => v != 0),6} | {100.0,5:F1} | {targetVals.Average(),8:E2} | {StdDev(targetVals),8:E2} | {targetVals.Min(),8:E2} | {targetVals.Max(),8:E2} | {"(target)",8} |        |");

            W("+---------------------------+--------+--------+-------+----------+----------+----------+----------+----------+--------+");
            W();
            W("Legend: SPARSE = <20% non-zero, SIGNAL = |r| > 0.15");

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

            var byGroup = _passing.GroupBy(r => (S(r, "ScanNumber"), S(r, "SourceFile"))).ToList();
            W($"Unique scans: {byGroup.Count}");
            W();

            var anomalous = new List<(string scan, string file, double sumTic, double maxTic, int count)>();

            foreach (var g in byGroup)
            {
                double sumTic = g.Sum(r => D(r, "TicNormalizedIntensity"));
                double maxTic = g.Max(r => D(r, "TicNormalizedIntensity"));

                if (sumTic > 0.5 || maxTic > 0.1)
                    anomalous.Add((g.Key.Item1, g.Key.Item2, sumTic, maxTic, g.Count()));
            }

            if (anomalous.Count > 0)
            {
                W($"Anomalous scans (SumTicNI > 0.5 or MaxTicNI > 0.1): {anomalous.Count}");
                W("+----------+------------------+----------+----------+-------+--------+");
                W("| Scan     | File             | SumTicNI | MaxTicNI | Count | Flag   |");
                W("+----------+------------------+----------+----------+-------+--------+");

                foreach (var (scan, file, sumTic, maxTic, count) in anomalous.OrderByDescending(x => x.sumTic).Take(20))
                {
                    string f = file.Length > 16 ? file.Substring(0, 13) + "..." : file;
                    string flag = sumTic > 0.5 ? "SUMHI" : "MAXHI";
                    W($"| {scan,-8} | {f,-16} | {sumTic,8:E2} | {maxTic,8:E2} | {count,5} | {flag,-6} |");
                }
                W("+----------+------------------+----------+----------+-------+--------+");
            }
            else
            {
                W("No anomalous scans detected. TIC normalization appears consistent.");
            }

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

            var ambiguous = _passing.Where(r => B(r, "IsIsobaricAmbiguous")).ToList();
            var unambiguous = _passing.Where(r => !B(r, "IsIsobaricAmbiguous")).ToList();

            W($"IsIsobaricAmbiguous = TRUE:  {ambiguous.Count:N0}");
            W($"IsIsobaricAmbiguous = FALSE: {unambiguous.Count:N0}");
            W();

            if (ambiguous.Count > 0 && unambiguous.Count > 0)
            {
                var ambTic = ambiguous.Select(r => D(r, "TicNormalizedIntensity")).ToList();
                var unambTic = unambiguous.Select(r => D(r, "TicNormalizedIntensity")).ToList();

                W("Comparison:");
                W($"  Ambiguous:   MeanTicNI={ambTic.Average():E4}, Median={Median(ambTic):E4}, MeanFragLen={ambiguous.Average(r => I(r, "FragmentLength")):F1}, FracBoth={100.0 * ambiguous.Count(r => B(r, "HasBothTerminalIons")) / ambiguous.Count:F1}%");
                W($"  Unambiguous: MeanTicNI={unambTic.Average():E4}, Median={Median(unambTic):E4}, MeanFragLen={unambiguous.Average(r => I(r, "FragmentLength")):F1}, FracBoth={100.0 * unambiguous.Count(r => B(r, "HasBothTerminalIons")) / unambiguous.Count:F1}%");
                W();

                if (ambTic.Average() > unambTic.Average() * 1.1)
                    W("NOTE: Ambiguous ions have higher mean intensity - may indicate real signals with mass degeneracy.");
                else if (unambTic.Average() > ambTic.Average() * 1.1)
                    W("NOTE: Unambiguous ions have higher mean intensity - ambiguous ions may be noisier.");
                else
                    W("NOTE: No systematic intensity difference between ambiguous and unambiguous ions.");
            }

            W();
            _out.Flush();
        }
    }
}