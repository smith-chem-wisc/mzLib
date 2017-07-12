using Chemistry;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net;
using UsefulProteomicsDatabases;

namespace Benchmark
{
    internal class Program
    {

        #region Private Methods

        private static void BenchmarkFormula()
        {
            Console.WriteLine("Starting benchmark BenchmarkFormula");

            int numRepetitions = 100000;

            Stopwatch stopWatch = new Stopwatch();

            var a = ChemicalFormula.ParseFormula("H1H{1}10 H{2}10 O20 O{16}20 O{17}20 O{18}20 C{12}100 C100 C{13}100 S{32}200 S200 S{33}200 S{34}200 S{36}200");
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var b = a.Formula + i;
            }
            stopWatch.Stop();
            Console.WriteLine("Time for getting formulas: " + stopWatch.Elapsed);

            Console.WriteLine("Benchmark BenchmarkFormula finished");
        }

        private static void BenchmarkFormula2()
        {
            Console.WriteLine("Starting benchmark BenchmarkFormula2");

            int numRepetitions = 100000;

            Stopwatch stopWatch = new Stopwatch();

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = ChemicalFormula.ParseFormula("H" + i + "H{1}10 H{2}10 O20 O{16}20 O{17}20 O{18}20 C{12}100 C100 C{13}100 S{32}200 S200 S{33}200 S{34}200 S{36}200");
                var b = a.Formula + i;
            }
            stopWatch.Stop();
            Console.WriteLine("Time for creating and getting formulas: " + stopWatch.Elapsed);

            Console.WriteLine("Benchmark BenchmarkFormula2 finished");
        }

        private static void BenchmarkGettingIsotopes()
        {
            Console.WriteLine("Starting benchmark BenchmarkGettingIsotopes");

            int numRepetitions = 10000000;

            Stopwatch stopWatch = new Stopwatch();

            long a = 0;
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                a += PeriodicTable.GetElement(20).Isotopes.Count();
            }
            stopWatch.Stop();
            Console.WriteLine("Time for getting isotopes1: " + stopWatch.Elapsed + " a = " + a);

            Console.WriteLine("Benchmark BenchmarkGettingIsotopes finished");
        }

        private static void BenchmarkIsotopicDistribution()
        {
            Console.WriteLine("Starting benchmark BenchmarkIsotopicDistribution");

            int numRepetitions = 100;

            Stopwatch stopWatch = new Stopwatch();

            var a = ChemicalFormula.ParseFormula("H100C100N100O100S100");
            double b = 0;
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                b += IsotopicDistribution.GetDistribution(a).Intensities.First();
            }
            stopWatch.Stop();
            Console.WriteLine("Time for generating isotopic distributions: " + stopWatch.Elapsed + " a = " + a);

            Console.WriteLine("Benchmark BenchmarkIsotopicDistribution finished");
        }

        private static void BenchmarkTimeGettingElementFromPeriodicTable()
        {
            Console.WriteLine("Starting benchmark BenchmarkTimeGettingElementFromPeriodicTable");

            int numRepetitions = 100000000;

            Stopwatch stopWatch = new Stopwatch();

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = PeriodicTable.GetElement(1);
                var b = a.Protons + a.AverageMass + 4;
            }
            stopWatch.Stop();
            Console.WriteLine("Time for getting by atomic number: " + stopWatch.Elapsed);

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = PeriodicTable.GetElement("H");
                var b = a.Protons + a.AverageMass + 4;
            }
            stopWatch.Stop();
            Console.WriteLine("Time for getting by atomic symbol: " + stopWatch.Elapsed);

            Console.WriteLine("Benchmark BenchmarkTimeGettingElementFromPeriodicTable finished");
        }

        private static void BenchmarkDatabaseLoadWrite()
        {
            Console.WriteLine("Starting benchmark BenchmarkDatabaseLoadWrite");

            Stopwatch stopWatch = new Stopwatch();

            Loaders.LoadElements("elements2.dat");
            IEnumerable<Modification> ya = PtmListLoader.ReadModsFromFile(@"ptmlist.txt").ToList();

            stopWatch.Restart();
            var a = ProteinDbLoader.LoadProteinXML(@"yeast_160126.xml.gz", true, ya, false, null, out Dictionary<string, Modification> um);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), a.Where(p => !p.IsDecoy).ToList(), "rewrite_yeast.xml");
            var b = ProteinDbLoader.LoadProteinXML(@"rewrite_yeast.xml", true, ya, false, null, out um);
            stopWatch.Stop();

            Console.WriteLine("Time for getting formulas: " + stopWatch.Elapsed);

            Console.WriteLine("Benchmark BenchmarkDatabaseLoadWrite finished");
        }

        private static void Main(string[] args)
        {
            Loaders.LoadElements("elements2.dat");

            int? minScan = null;
            int? maxScan = null;
            double deconvolutionTolerancePpm = 5;
            double aggregationTolerancePpm = 5;
            //string filename = @"C:\Users\stepa\Desktop\DeconvolutionStuff\03-01-17_B2A_targeted_td_yeast_fract6_intact.raw";
            //string filename = @"C:\Users\stepa\Desktop\neucode_for_deconv\10-23-15_A_fract5_rep2.raw";
            //string filename = @"C:\Users\stepa\Desktop\neucode_for_deconv\10-28-15_S_fract5_rep1.raw";
            string filename = @"C:\Users\stepa\Desktop\RobDeconvolution\B02_06_161103_A1_HCD_OT_4ul.raw";
            int maxAssumedChargeState = 10;
            double intensityRatioLimit = 4;
            Func<IMzPeak, bool> peakFilter = b => (b as ThermoMzPeak).SignalToNoise > 2;

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> f = ThermoStaticData.LoadAllStaticData(filename);
            //IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> f = ThermoDynamicData.InitiateDynamicConnection(filename);

            Func<IMsDataScan<IMzSpectrum<IMzPeak>>, bool> scanFilterFunc = b => b.MsnOrder == 1;
            IEnumerable<DeconvolutionFeatureWithMassesAndScans> nice = f.Deconvolute(minScan, maxScan, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit, peakFilter, aggregationTolerancePpm, scanFilterFunc);

            Console.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.NumPeaks).Take(10)));

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"out-numPeaks.txt"))
            {
                file.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.NumPeaks).Select(b => b.OneLineString())));
            }

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"out-int.txt"))
            {
                file.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.TotalIntensity).Select(b => b.OneLineString())));
            }

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"out-Mass.txt"))
            {
                file.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.Mass).Select(b => b.OneLineString())));
            }

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"out-MinScanIndex.txt"))
            {
                file.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => b.MinScanIndex).Select(b => b.OneLineString())));
            }



            using (WebClient Client = new WebClient())
                Client.DownloadFile(@"http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some", "Dddd.temp");

            DoubleRange r = new DoubleRange(-187, double.PositiveInfinity);
            Console.WriteLine(r);
            Console.WriteLine(r.ToString());

            //Dictionary<string, Modification> um;
            //ProteinDbLoader.LoadProteinXML(@"C:\Users\stepa\Desktop\01012017_MM_ONLYGPTMD.xml", true, new List<Modification>(), false, new List<string> { "GO", "EnsemblFungi" }, null, out um);

            //ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Desktop\PrecursorProblems\2016_080902_SMC_EC_Glyco_EThcD.raw");

            //Mzml.LoadAllStaticData(@"C:\Users\stepa\Source\Repos\MetaMorpheus\Test\bin\Debug\ok.mzML");
            //var oddk = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac1.raw");
            //ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Desktop\yeast_tmt\m04667.raw");
            //ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Desktop\human_spike\C14-11130.raw");
            //ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac18.raw");
            ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Desktop\CoIsolation\05-11-17_YL_25iso.raw");

            //var hheh = oddk.GetMsScansInTimeRange(47.2469, 47.25693).ToList();

            //Mzml.LoadAllStaticData(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac28.mzML");
            //var okff = ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac28.raw");
            //var okff = ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Data\ForRyan\golden.raw");
            //MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(okff, @"C:\Users\stepa\Data\ForRyan\adsfjk.mzML");

            //Mzml.LoadAllStaticData(@"C:\Users\stepa\Desktop\02-15-17_Cys-tag_light\02-14-17_Cl-1_rep1.mzML");
            //using (var nice = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Desktop\02-15-17_Cys-tag_light\02-14-17_Cl-1_rep1.raw"))
            //{
            //    Console.WriteLine(nice.GetOneBasedScan(1000).RetentionTime);
            //}

            //// OLD MASS SPEC
            //var theFiles = new List<string>{
            //    //@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac17.raw",
            //    //@"C:\Users\stepa\Data\CalibrationPaperData\Mouse\04-30-13_CAST_Frac5_4uL.raw",
            //    //@"C:\Users\stepa\Data\CalibrationPaperData\Yeast\12-10-16_A17A_yeast_BU_fract9_rep1_8uL.raw",
            //    //@"C:\Users\stepa\Desktop\MvsMM\04-21-17_Lys_1-200_rep1.raw",
            //    //@"C:\Users\stepa\Desktop\MvsMM\04-21-17_Lys_1-200_rep1.mzML",
            //    @"C:\Users\stepa\Data\CalibrationPaperData\Mouse\04-29-13_B6_Frac7_5uL.raw"
            //};

            //// Params
            //var tols = new List<Tolerance> { new Tolerance("5 PPM") };
            //var isotopeRatios = new List<int> { 4 };
            //var maxAssumedChargeState = 10;

            //foreach (var theFile in theFiles)
            //{
            //    var okff = ThermoStaticData.LoadAllStaticData(theFile);
            //    //var okff = Mzml.LoadAllStaticData(theFile);

            //    int countScans = 0;
            //    int totalHaveMMandCharge = 0;

            //    var totalHaveMyMass = new int[1, 1];
            //    var totalMatch = new int[1, 1];
            //    //foreach (var scanWithPrec in okff.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            //    var scanWithPrec = okff.GetOneBasedScan(11042) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            //    {
            //        countScans++;

            //        Console.WriteLine("Scan " + scanWithPrec.OneBasedScanNumber + " ; isolation=" + scanWithPrec.IsolationMz + " ; mm=" + scanWithPrec.SelectedIonMonoisotopicGuessMz + " ; charge=" + scanWithPrec.SelectedIonChargeStateGuess);

            //        if (scanWithPrec.SelectedIonMonoisotopicGuessMz.HasValue && scanWithPrec.SelectedIonChargeStateGuess.HasValue)
            //        {
            //            totalHaveMMandCharge++;
            //        }

            //        for (int i = 0; i < tols.Count; i++)
            //        {
            //            var tol = tols[i];
            //            for (int j = 0; j < isotopeRatios.Count; j++)
            //            {
            //                var isotopeRatio = isotopeRatios[j];
            //                var mzEnvelopesWithCharges = scanWithPrec.GetIsolatedMassesAndCharges(okff.GetOneBasedScan(scanWithPrec.OneBasedPrecursorScanNumber).MassSpectrum, maxAssumedChargeState, tol, isotopeRatio).ToList();

            //                if (mzEnvelopesWithCharges.Count() > 0)
            //                    totalHaveMyMass[i, j]++;

            //                if (scanWithPrec.SelectedIonMonoisotopicGuessMz.HasValue && scanWithPrec.SelectedIonChargeStateGuess.HasValue)
            //                {
            //                    if (mzEnvelopesWithCharges.Any(bd => tol.Within(bd.Item1.First().Mz.ToMass(bd.Item2), scanWithPrec.SelectedIonMonoisotopicGuessMz.Value.ToMass(scanWithPrec.SelectedIonChargeStateGuess.Value))))
            //                    {
            //                        totalMatch[i, j]++;
            //                        Console.WriteLine("Match!");
            //                    }
            //                    else
            //                    {
            //                        Console.WriteLine(string.Join(Environment.NewLine, mzEnvelopesWithCharges.Select(b => "\t" + b.Item2 + " : " + string.Join(",", b.Item1))));
            //                    }
            //                }
            //                else
            //                {
            //                    Console.WriteLine(string.Join(Environment.NewLine, mzEnvelopesWithCharges.Select(b => "\t" + b.Item2 + " : " + string.Join(",", b.Item1))));
            //                }
            //            }
            //        }
            //    }

            //    Console.WriteLine("countScans: " + countScans);
            //    Console.WriteLine("totalHaveMMandCharge: " + totalHaveMMandCharge);

            //    for (int i = 0; i < tols.Count; i++)
            //    {
            //        var tol = tols[i];
            //        for (int j = 0; j < isotopeRatios.Count; j++)
            //        {
            //            Console.WriteLine("i = " + i + " j = " + j);
            //            Console.WriteLine("totalHaveMyMass: " + totalHaveMyMass[i, j]);
            //            Console.WriteLine("totalMatch: " + totalMatch[i, j]);
            //        }
            //    }
            //}

            //using (var nice = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Data\CalibrationPaperData\Mouse\04-30-13_CAST_Frac5_4uL.raw"))
            ////{
            //using (var nice = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac17.raw"))
            //{
            //var hehdfe = nice.GetOneBasedScan(168) as IMsDataScanWithPrecursor<ThermoSpectrum>;
            //Console.WriteLine("Scan " + hehdfe.OneBasedScanNumber + " ; isolation=" + hehdfe.IsolationMz + " ; mm=" + hehdfe.SelectedIonMonoisotopicGuessMz);

            //var fdf = hehdfe.GetIsolatedMassesAndCharges(nice.GetOneBasedScan(hehdfe.OneBasedPrecursorScanNumber).MassSpectrum, 10, new Tolerance("20 PPM"), 10, 1);

            //Console.WriteLine(fdf.Count() + ";" + hehdfe.SelectedIonMonoisotopicGuessMz.HasValue);

            //Console.WriteLine(string.Join(Environment.NewLine, fdf.Select(b => "\t" + b.Item1 + "; " + b.Item2 + "; " + b.Item1.First())));

            //    Console.WriteLine();

            //    int totalHaveMM = 0;
            //    int totalHaveMMandMatch = 0;
            //    int totalHaveMyMass = 0;
            //    int totalHaveMMandMatchAll = 0;
            //    int totalHaveMyMassAll = 0;
            //    foreach (var hehdfe in nice.OfType<IMsDataScanWithPrecursor<ThermoSpectrum>>())
            //    {
            //        Console.WriteLine("Scan " + hehdfe.OneBasedScanNumber + " ; isolation=" + hehdfe.IsolationMz + " ; mm=" + hehdfe.SelectedIonMonoisotopicGuessMz + " ; charge=" + hehdfe.SelectedIonChargeStateGuess);

            //        //var fdf = hehdfe.GetIsolatedMassesAndCharges(nice.GetOneBasedScan(hehdfe.OneBasedPrecursorScanNumber).MassSpectrum, 10, tol, 10, 1).ToList();

            //        //if (fdf.Count() > 0)
            //        //    totalHaveMyMass++;

            //        //if (hehdfe.SelectedIonMonoisotopicGuessMz.HasValue)
            //        //{
            //        //    totalHaveMM++;
            //        //    if (fdf.Any(bd => tol.Within(bd.Item1.First().ToMass(bd.Item2), hehdfe.SelectedIonMonoisotopicGuessMz.Value.ToMass(hehdfe.SelectedIonChargeStateGuess.Value))))
            //        //        totalHaveMMandMatch++;
            //        //}

            //        if (fdfAll.Count() > 0)
            //            totalHaveMyMassAll++;
            //        Console.WriteLine(fdfAll.Count() + ";" + hehdfe.SelectedIonMonoisotopicGuessMz.HasValue);

            //        if (hehdfe.SelectedIonMonoisotopicGuessMz.HasValue)
            //        {
            //            totalHaveMM++;
            //            var massFromScan = hehdfe.SelectedIonMonoisotopicGuessMz.Value.ToMass(hehdfe.SelectedIonChargeStateGuess.Value);
            //            if (fdfAll.Any(bd => tol.Within(bd.Item1.First().ToMass(bd.Item2), massFromScan)))
            //                totalHaveMMandMatchAll++;

            //            if (totalHaveMMandMatchAll - totalHaveMyMassAll < -6)
            //                    Console.WriteLine(totalHaveMMandMatchAll + " not equal " + totalHaveMyMassAll);
            //        }

            //        Console.WriteLine();
            //    }
            //}

            //using (var nice = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Data\CalibrationPaperData\Jurkat\120426_Jurkat_highLC_Frac17.raw"))
            //{
            //    var ok = nice.GetOneBasedScan(24676);
            //    var hm = ok as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            //    hm.RecomputeSelectedPeak(nice.GetOneBasedScan(hm.OneBasedPrecursorScanNumber).MassSpectrum);
            //}

            //using (var nice = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Desktop\02-15-17_Cys-tag_light\02-14-17_Cl-1_rep1.raw"))
            //{
            //    var ok = nice.GetOneBasedScan(71291);
            //    var hm = ok as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;

            //    var prevSpectrum = nice.GetOneBasedScan(hm.OneBasedPrecursorScanNumber).MassSpectrum;

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.RecomputeChargeState(prevSpectrum, 0.01, 4);

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.ComputeSelectedPeakIntensity(prevSpectrum);

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.ComputeMonoisotopicPeakIntensity(prevSpectrum);

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.RecomputeSelectedPeak(prevSpectrum);

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.RecomputeMonoisotopicPeak(prevSpectrum, 0.01, 0.3);

            //    Console.WriteLine(hm.SelectedIonGuessChargeStateGuess + Environment.NewLine
            //        + hm.SelectedIonGuessMZ + Environment.NewLine
            //        + hm.SelectedIonGuessIntensity + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicMZ + Environment.NewLine
            //        + hm.SelectedIonGuessMonoisotopicIntensity + Environment.NewLine);

            //    hm.RecomputeSelectedPeak(nice.GetOneBasedScan(hm.OneBasedPrecursorScanNumber).MassSpectrum);
            //}

            PopulatePeriodicTable();

            BenchmarkFormula();
            Console.WriteLine("");
            BenchmarkFormula2();
            Console.WriteLine("");
            BenchmarkTimeGettingElementFromPeriodicTable();
            Console.WriteLine("");
            BenchmarkGettingIsotopes();
            Console.WriteLine("");
            BenchmarkIsotopicDistribution();
            Loaders.LoadElements(@"elements.tmp");
            Console.WriteLine("");
            BenchmarkDatabaseLoadWrite();
        }

        private static void PopulatePeriodicTable()
        {
            var elementH = new Element("H", 1, 1.007975);
            PeriodicTable.Add(elementH);
            elementH.AddIsotope(1, 1.00782503223, 0.999885);
            elementH.AddIsotope(2, 2.01410177812, 0.000115);

            var elementC = new Element("C", 6, 12.0106);
            PeriodicTable.Add(elementC);
            elementC.AddIsotope(12, 12, 0.9893);
            elementC.AddIsotope(13, 13.00335483507, 0.0107);

            var elementN = new Element("N", 7, 14.006855);
            PeriodicTable.Add(elementN);
            elementN.AddIsotope(15, 15.00010889888, 0.00364);
            elementN.AddIsotope(14, 14.00307400443, 0.99636);

            var elementO = new Element("O", 8, 15.9994);
            PeriodicTable.Add(elementO);
            elementO.AddIsotope(16, 15.99491461957, 0.99757);
            elementO.AddIsotope(17, 16.99913175650, 0.00038);
            elementO.AddIsotope(18, 17.99915961286, 0.00205);

            var elementFe = new Element("Fe", 26, 55.845);
            PeriodicTable.Add(elementFe);
            elementFe.AddIsotope(54, 53.93960899, 0.05845);
            elementFe.AddIsotope(56, 55.93493633, 0.91754);
            elementFe.AddIsotope(57, 56.93539284, 0.02119);
            elementFe.AddIsotope(58, 57.93327443, 0.00282);

            var elementBr = new Element("Br", 35, 79.904);
            PeriodicTable.Add(elementBr);
            elementBr.AddIsotope(79, 78.9183376, 0.5069);
            elementBr.AddIsotope(81, 80.9162897, 0.4931);

            var elementCa = new Element("Ca", 20, 40.078);
            PeriodicTable.Add(elementCa);
            elementCa.AddIsotope(40, 39.962590863, 0.96941);
            elementCa.AddIsotope(42, 41.95861783, 0.00647);
            elementCa.AddIsotope(43, 42.95876644, 0.00135);
            elementCa.AddIsotope(44, 43.95548156, 0.02086);
            elementCa.AddIsotope(46, 45.9536890, 0.00004);
            elementCa.AddIsotope(48, 47.95252276, 0.00187);

            var elementS = new Element("S", 16, 32.0675);
            PeriodicTable.Add(elementS);
            elementS.AddIsotope(32, 31.9720711744, 0.9499);
            elementS.AddIsotope(33, 32.9714589098, 0.0075);
            elementS.AddIsotope(34, 33.967867004, 0.0425);
            elementS.AddIsotope(36, 35.96708071, 0.0001);

            var elementSe = new Element("Se", 34, 78.971);
            PeriodicTable.Add(elementSe);
            elementSe.AddIsotope(74, 73.922475934, 0.0089);
            elementSe.AddIsotope(76, 75.919213704, 0.0937);
            elementSe.AddIsotope(77, 76.919914154, 0.0763);
            elementSe.AddIsotope(78, 77.91730928, 0.2377);
            elementSe.AddIsotope(80, 79.9165218, 0.4961);
            elementSe.AddIsotope(82, 81.9166995, 0.0873);

            var elementZr = new Element("Zr", 40, 91.224);
            PeriodicTable.Add(elementZr);
            elementZr.AddIsotope(90, 89.9046977, 0.5145);
            elementZr.AddIsotope(91, 90.9056396, 0.1122);
            elementZr.AddIsotope(92, 91.9050347, 0.1715);
            elementZr.AddIsotope(94, 93.9063108, 0.1738);
        }

        #endregion Private Methods

    }
}