using Chemistry;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace ConsoleApp1
{
    internal class Program
    {
        #region Private Methods

        private static void Main(string[] args)
        {
            Loaders.LoadElements("elements.dat");

            var sameMassTolerancePpm = new PpmTolerance(10);

            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"C:\Users\stepa\Desktop\Decon\04-29-13_B6_Frac5_4uL.raw");

            using (StreamWriter output = new StreamWriter(@"04-29-13_B6_Frac5_4uL.tsv"))
            {
                var intensityRatio = 5;
                var deconvolutionTolerancePpm = 20;
                var maxAssumedChargeState = 10;
                Console.WriteLine("intensityRatio, deconvolutionTolerancePpm, maxAssumedChargeState: " + (intensityRatio, deconvolutionTolerancePpm, maxAssumedChargeState));
                Tolerance deconvolutionTolerance = new PpmTolerance(deconvolutionTolerancePpm);

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

                output.WriteLine("Scan\tType");
                //foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().Where(b => b.OneBasedScanNumber == 3356))
                foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>().Where(b => b.OneBasedScanNumber >= 3514))
                //foreach (var ok in a.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
                {
                    double thermoMz = double.NaN;
                    if (ok.SelectedIonChargeStateGuess.HasValue)
                    {
                        thermoMz = ok.SelectedIonMZ;
                        if (ok.SelectedIonMonoisotopicGuessMz.HasValue)
                            thermoMz = ok.SelectedIonMonoisotopicGuessMz.Value;
                        output.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + thermoMz.ToString("G6"));
                        output.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + thermoMz.ToMass(ok.SelectedIonChargeStateGuess.Value).ToString("G6"));
                    }
                    else
                    {
                        output.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + " ");
                        output.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + " ");
                    }

                    var precursorSpectrum = a.GetOneBasedScan(ok.OneBasedPrecursorScanNumber.Value).MassSpectrum;

                    var oldResults = ok.GetIsolatedMassesAndChargesOld(precursorSpectrum, maxAssumedChargeState, deconvolutionTolerance, intensityRatio).OrderBy(b => b.Item1.First().Mz.ToMass(b.Item2)).ToList();

                    output.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMzs\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToString("G6"))));
                    output.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMasses\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToMass(b.Item2).ToString("G6"))));

                    var newResults = ok.GetIsolatedMassesAndCharges(precursorSpectrum, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).OrderBy(b => b.monoisotopicMass).ToList();

                    output.WriteLine(ok.OneBasedScanNumber + "\tNewLowestObservedMzs\t" + string.Join("\t", newResults.Select(b => b.peaks.OrderBy(c => c.Item1).First().Item1.ToString("G6"))));
                    output.WriteLine(ok.OneBasedScanNumber + "\tNewMonoisotopicMasses\t" + string.Join("\t", newResults.Select(b => b.monoisotopicMass.ToString("G6"))));

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
                        {
                            numMatchToOnlyOld++;
                        }
                        if (matchToNew && !matchToOld)
                            numMatchToOnlyNew++;
                        if (matchToOld && matchToNew)
                            numMatchToBoth++;
                        if (!matchToOld && !matchToNew && thermoMz > 500)
                        {
                            Console.WriteLine(ok.OneBasedScanNumber + "\tThermoMz\t" + thermoMz.ToString("G6"));
                            Console.WriteLine(ok.OneBasedScanNumber + "\tThermoMass\t" + thermoMz.ToMass(ok.SelectedIonChargeStateGuess.Value).ToString("G6"));
                            Console.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMzs\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToString("G6"))));
                            Console.WriteLine(ok.OneBasedScanNumber + "\tOldMonoisotopicMasses\t" + string.Join("\t", oldResults.Select(b => b.Item1.First().Mz.ToMass(b.Item2).ToString("G6"))));
                            Console.WriteLine(ok.OneBasedScanNumber + "\tNewLowestObservedMzs\t" + string.Join("\t", newResults.Select(b => b.peaks.OrderBy(c => c.Item1).First().Item1.ToString("G6"))));
                            Console.WriteLine(ok.OneBasedScanNumber + "\tNewMonoisotopicMasses\t" + string.Join("\t", newResults.Select(b => b.monoisotopicMass.ToString("G6"))));
                        }
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

        #endregion Private Methods
    }
}