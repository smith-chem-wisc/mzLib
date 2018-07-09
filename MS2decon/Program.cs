using Fclp;
using IO.MzML;
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

namespace MS2decon
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            Loaders.LoadElements("elements.dat");

            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.DeconvolutionTolerancePpm)
             .As("DeconvolutionTolerancePpm");

            p.Setup(arg => arg.MinScan)
             .As("MinScan");

            p.Setup(arg => arg.MaxScan)
             .As("MaxScan");

            p.Setup(arg => arg.MinAssumedChargeState)
             .As("MinAssumedChargeState");

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
                MsDataFile myMsDataFile;
                if (Path.GetExtension(p.Object.FilePath).Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                    myMsDataFile = Mzml.LoadAllStaticData(p.Object.FilePath);
                else
                    myMsDataFile = ThermoStaticData.LoadAllStaticData(p.Object.FilePath);

                using (StreamWriter output = new StreamWriter(@"MS2DeconvolutionOutput-" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture) + ".tsv"))
                {
                    output.WriteLine("Mass\tNumPeaks\tNumScans\tMinScan\tMaxScan\tAverageElutionTime\tIntensity\tObservedCharges\tMostIntenseCharge\tMostIntenseMz\tNumPeaksInMostIntenseEnvelope");
                    foreach (var ok in myMsDataFile.GetAllScansList().Where(x => x.MsnOrder != 1))
                    {
                        if ((!p.Object.MinScan.HasValue || ok.OneBasedScanNumber >= p.Object.MinScan) && (!p.Object.MaxScan.HasValue || ok.OneBasedScanNumber <= p.Object.MaxScan))
                        {
                            var hmm = ok.MassSpectrum.Deconvolute(new MzRange(0, double.PositiveInfinity), p.Object.MinAssumedChargeState, p.Object.MaxAssumedChargeState, p.Object.DeconvolutionTolerancePpm, p.Object.IntensityRatioLimit).ToList();

                            List<DeconvolutionFeatureWithMassesAndScans> currentListOfGroups = new List<DeconvolutionFeatureWithMassesAndScans>();

                            foreach (var isotopicEnvelope in hmm)
                            {
                                DeconvolutionFeatureWithMassesAndScans matchingGroup = null;
                                var mass = isotopicEnvelope.monoisotopicMass;
                                foreach (var possibleGroup in currentListOfGroups)
                                {
                                    var possibleGroupMass = possibleGroup.Mass;
                                    if (Math.Abs(mass - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass + 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass + 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass + 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass - 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass - 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm ||
                                        Math.Abs(mass - 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= p.Object.AggregationTolerancePpm)
                                    {
                                        matchingGroup = possibleGroup;
                                        matchingGroup.AddEnvelope(isotopicEnvelope, ok.OneBasedScanNumber, ok.RetentionTime);
                                        break;
                                    }
                                }

                                if (matchingGroup == null)
                                {
                                    var newGroupScans = new DeconvolutionFeatureWithMassesAndScans();
                                    newGroupScans.AddEnvelope(isotopicEnvelope, ok.OneBasedScanNumber, ok.RetentionTime);
                                    currentListOfGroups.Add(newGroupScans);
                                }
                            }
                            foreach (var ook in currentListOfGroups)
                                output.WriteLine(ook.OneLineString());
                        }
                    }
                }
            }
            else
            {
                Console.WriteLine("BAD PARAMETERS");
                Console.WriteLine(result.ErrorText);
            }
        }
    }

    internal class ApplicationArguments
    {
        public int? MinScan { get; set; } = null;
        public int? MaxScan { get; set; } = null;
        public int MinAssumedChargeState { get; set; } = 1;
        public int MaxAssumedChargeState { get; set; } = 10;
        public double DeconvolutionTolerancePpm { get; set; } = 20;
        public double IntensityRatioLimit { get; set; } = 5;
        public double AggregationTolerancePpm { get; set; } = 5;
        public string FilePath { get; set; }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("FilePath: " + FilePath);
            sb.AppendLine("MinScan: " + MinScan);
            sb.AppendLine("MaxScan: " + MaxScan);
            sb.AppendLine("MinAssumedChargeState: " + MinAssumedChargeState);
            sb.AppendLine("MaxAssumedChargeState: " + MaxAssumedChargeState);
            sb.AppendLine("DeconvolutionTolerancePpm: " + DeconvolutionTolerancePpm);
            sb.AppendLine("IntensityRatioLimit: " + IntensityRatioLimit);
            sb.AppendLine("AggregationTolerancePpm: " + AggregationTolerancePpm);
            return sb.ToString();
        }
    }
}