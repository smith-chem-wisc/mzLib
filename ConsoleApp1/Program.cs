using Chemistry;
using Fclp;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace ConsoleApp1
{
    internal class Program
    {
        #region Private Methods

        private static void Main(string[] args)
        {
            Loaders.LoadElements("elements.dat");

            if (args.Length > 0)
            {
                DoFileDecon(args);
            }
            else { DeconTest(); }
        }

        private static void DeconTest()
        {
            List<string> files = new List<string>
            {
                @"C:\Users\stepa\Desktop\Decon\04-29-13_B6_Frac5_4uL.raw",
                @"C:\Users\stepa\Desktop\Decon\12-10-16_A17A_yeast_BU_fract9_rep1_8uL.raw",
                @"C:\Users\stepa\Desktop\Decon\120426_Jurkat_highLC_Frac17.raw",
                @"C:\Users\stepa\Desktop\Decon\07-26-17_YL_150mins_steeper-grad_yeast.raw",
                @"C:\Users\stepa\Desktop\Decon\07-28-17_19_rep2_mouse.raw",
                @"C:\Users\stepa\Desktop\Decon\09-01-17_iodo_1-4th_rep1_human.raw",
            };

            var sameMassTolerancePpm = new PpmTolerance(10);

            var intensityRatio = 5;
            var deconvolutionTolerancePpm = 20;
            var maxAssumedChargeState = 10;

            Tolerance deconvolutionTolerance = new PpmTolerance(deconvolutionTolerancePpm);

            foreach (var file in files)
            {
                Console.WriteLine(file);
                ThermoStaticData a = ThermoStaticData.LoadAllStaticData(file);
                //ThermoDynamicData a = ThermoDynamicData.InitiateDynamicConnection(@"C:\Users\stepa\Desktop\Decon\09-01-17_iodo_1-4th_rep1_human.raw");

                //using (StreamWriter output = new StreamWriter(@"output.tsv"))
                //{
                //Console.WriteLine("intensityRatio, deconvolutionTolerancePpm, maxAssumedChargeState: " + (intensityRatio, deconvolutionTolerancePpm, maxAssumedChargeState));

                //int goodScans = 0;

                int numMatchToOnlyOld = 0;
                int numMatchToOnlyNew = 0;
                int numMatchToBoth = 0;
                int numScansWithThermoMasses = 0;

                int numMatchToOld = 0;
                int numMatchToOldPlusMM = 0;
                int numMatchToOldMinusMM = 0;
                int numMatchToNew = 0;
                int numMatchToNewPlusMM = 0;
                int numMatchToNewMinusMM = 0;

                //output.WriteLine("Scan\tType");
                //foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().Where(b => b.OneBasedScanNumber >= 8356))
                //foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().Skip(2000))
                foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
                {
                    double thermoMz = double.NaN;
                    if (ok.SelectedIonChargeStateGuess.HasValue)
                    {
                        thermoMz = ok.SelectedIonMZ;
                        if (ok.SelectedIonMonoisotopicGuessMz.HasValue)
                            thermoMz = ok.SelectedIonMonoisotopicGuessMz.Value;
                        //output.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + thermoMz.ToString("G6"));
                        //output.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + thermoMz.ToMass(ok.SelectedIonChargeStateGuess.Value).ToString("G6"));
                    }
                    else
                    {
                        //output.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + " ");
                        //output.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + " ");
                    }

                    var precursorSpectrum = a.GetOneBasedScan(ok.OneBasedPrecursorScanNumber.Value).MassSpectrum;

                    var oldResults = ok.GetIsolatedMassesAndChargesOld(precursorSpectrum, maxAssumedChargeState, deconvolutionTolerance, intensityRatio).OrderBy(b => b.Item1.First().Mz.ToMass(b.Item2)).ToList();

                    //output.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMzs\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToString("G6"))));
                    //output.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMasses\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToMass(b.Item2).ToString("G6"))));

                    var newResults = ok.GetIsolatedMassesAndCharges(precursorSpectrum, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).OrderBy(b => b.monoisotopicMass).ToList();

                    //output.WriteLine(ok.OneBasedScanNumber + "\tNewLowestObservedMzs\t" + string.Join("\t", newResults.Select(b => b.peaks.OrderBy(c => c.Item1).First().Item1.ToString("G6"))));
                    //output.WriteLine(ok.OneBasedScanNumber + "\tNewMonoisotopicMasses\t" + string.Join("\t", newResults.Select(b => b.monoisotopicMass.ToString("G6"))));

                    if (ok.SelectedIonChargeStateGuess.HasValue)
                    {
                        numScansWithThermoMasses++;
                        var thermoMass = thermoMz.ToMass(ok.SelectedIonChargeStateGuess.Value);
                        bool matchToOld = false;
                        bool matchToNew = false;
                        bool matchToOldPlusMM = false;
                        bool matchToNewPlusMM = false;
                        bool matchToOldMinusMM = false;
                        bool matchToNewMinusMM = false;
                        foreach (var hm in oldResults)
                        {
                            if (sameMassTolerancePpm.Within(hm.Item1.First().Mz.ToMass(hm.Item2), thermoMass))
                            {
                                matchToOld = true;
                            }
                            if (sameMassTolerancePpm.Within(hm.Item1.First().Mz.ToMass(hm.Item2) + 1.0029, thermoMass))
                            {
                                matchToOldPlusMM = true;
                            }
                            if (sameMassTolerancePpm.Within(hm.Item1.First().Mz.ToMass(hm.Item2) - 1.0029, thermoMass))
                            {
                                matchToOldMinusMM = true;
                            }
                        }
                        foreach (var hm in newResults)
                        {
                            if (sameMassTolerancePpm.Within(hm.monoisotopicMass, thermoMass))
                            {
                                matchToNew = true;
                            }
                            if (sameMassTolerancePpm.Within(hm.monoisotopicMass + 1.0029, thermoMass))
                            {
                                matchToNewPlusMM = true;
                            }
                            if (sameMassTolerancePpm.Within(hm.monoisotopicMass - 1.0029, thermoMass))
                            {
                                matchToNewMinusMM = true;
                            }
                        }
                        numMatchToOld += matchToOld ? 1 : 0;
                        numMatchToOldPlusMM += matchToOldPlusMM ? 1 : 0;
                        numMatchToOldMinusMM += matchToOldMinusMM ? 1 : 0;
                        numMatchToNew += matchToNew ? 1 : 0;
                        numMatchToNewPlusMM += matchToNewPlusMM ? 1 : 0;
                        numMatchToNewMinusMM += matchToNewMinusMM ? 1 : 0;
                        if (matchToOld && !matchToNew)
                            numMatchToOnlyOld++;
                        if (matchToNew && !matchToOld)
                            numMatchToOnlyNew++;
                        if (matchToOld && matchToNew)
                            numMatchToBoth++;
                        //if (matchToOld && !matchToNew && thermoMz > 500)
                        //{
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + thermoMz.ToString("G6"));
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + thermoMz.ToMass(ok.SelectedIonChargeStateGuess.Value).ToString("G6"));
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMzs\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToString("G6"))));
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMasses\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToMass(b.Item2).ToString("G6"))));
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tNewLowestObservedMzs\t" + string.Join("\t", newResults.Select(b => b.peaks.OrderBy(c => c.Item1).First().Item1.ToString("G6"))));
                        //    Console.WriteLine(ok.OneBasedScanNumber + "\tNewMonoisotopicMasses\t" + string.Join("\t", newResults.Select(b => b.monoisotopicMass.ToString("G6"))));

                        //}
                        //else
                        //{
                        //    Console.WriteLine(ok.OneBasedScanNumber + " good ");
                        //}
                    }
                }

                Console.WriteLine("numMatchToOld: " + numMatchToOld);
                Console.WriteLine("numMatchToOldPlusMM: " + numMatchToOldPlusMM);
                Console.WriteLine("numMatchToOldMinusMM: " + numMatchToOldMinusMM);
                Console.WriteLine("numMatchToNew: " + numMatchToNew);
                Console.WriteLine("numMatchToNewPlusMM: " + numMatchToNewPlusMM);
                Console.WriteLine("numMatchToNewMinusMM: " + numMatchToNewMinusMM);

                Console.WriteLine("numMatchToOnlyOld: " + numMatchToOnlyOld);
                Console.WriteLine("numMatchToOnlyNew: " + numMatchToOnlyNew);
                Console.WriteLine("numMatchToBoth: " + numMatchToBoth);
                Console.WriteLine("numScansWithThermoMasses: " + numScansWithThermoMasses);
            }
        }

        private static void DoFileDecon(string[] args)
        {
            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.AggregationTolerancePpm)
             .As("AggregationTolerancePpm");

            p.Setup(arg => arg.DeconvolutionTolerancePpm)
             .As("DeconvolutionTolerancePpm");

            p.Setup(arg => arg.MinScan)
             .As("MinScan");

            p.Setup(arg => arg.MaxScan)
             .As("MaxScan");

            p.Setup(arg => arg.MaxAssumedChargeState)
             .As("MaxAssumedChargeState");

            p.Setup(arg => arg.IntensityRatioLimit)
             .As("IntensityRatioLimit");

            p.Setup(arg => arg.FilePath)
             .As("FilePath").
             Required();

            var result = p.Parse(args);

            Console.WriteLine("Running deconvolution using the following parameters:");
            Console.WriteLine(p.Object);

            if (result.HasErrors == false)
            {
                ThermoStaticData a = ThermoStaticData.LoadAllStaticData(p.Object.FilePath);

                using (StreamWriter output = new StreamWriter(@"DeconvolutionOutput-" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture) + ".tsv"))
                {
                    output.WriteLine("Mass\tNumPeaks\tNumScans\tMinScan\tMaxScan\tAverageElutionTime\tIntensity\tMostIntenseCharge\tMostIntenseMz\tNumPeaksInMostIntenseEnvelope");
                    foreach (var nice in a.Deconvolute(p.Object.MinScan, p.Object.MaxScan, p.Object.MaxAssumedChargeState, p.Object.DeconvolutionTolerancePpm, p.Object.IntensityRatioLimit, p.Object.AggregationTolerancePpm, b => b.MsnOrder == 1).OrderByDescending(b => b.TotalIntensity))
                    {
                        output.WriteLine(nice.OneLineString());
                    }
                }
            }
            else
            {
                Console.WriteLine("BAD PARAMETERS");
                Console.WriteLine(result.ErrorText);
            }
        }

        #endregion Private Methods
    }

    internal class ApplicationArguments
    {
        #region Public Properties

        public int? MinScan { get; set; } = null;
        public int? MaxScan { get; set; } = null;
        public int MaxAssumedChargeState { get; set; } = 10;
        public double DeconvolutionTolerancePpm { get; set; } = 20;
        public double IntensityRatioLimit { get; set; } = 5;
        public double AggregationTolerancePpm { get; set; } = 5;
        public string FilePath { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("FilePath: " + FilePath);
            sb.AppendLine("MinScan: " + MinScan);
            sb.AppendLine("MaxScan: " + MaxScan);
            sb.AppendLine("MaxAssumedChargeState: " + MaxAssumedChargeState);
            sb.AppendLine("DeconvolutionTolerancePpm: " + DeconvolutionTolerancePpm);
            sb.AppendLine("IntensityRatioLimit: " + IntensityRatioLimit);
            sb.AppendLine("AggregationTolerancePpm: " + AggregationTolerancePpm);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}