using Fclp;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace ConsoleApp1
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            Loaders.LoadElements("elements.dat");

            DoFileDecon(args);
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

            p.Setup(arg => arg.MinAssumedChargeState)
             .As("MinAssumedChargeState");

            p.Setup(arg => arg.MaxAssumedChargeState)
             .As("MaxAssumedChargeState");

            p.Setup(arg => arg.IntensityRatioLimit)
             .As("IntensityRatioLimit");

            p.Setup(arg => arg.AverageScans)
             .As("AverageScans");

            p.Setup(arg => arg.NumScansRequired)
             .As("NumScansRequired");

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

                if (p.Object.AverageScans > 1)
                {
                    myMsDataFile = new SummedMsDataFile(myMsDataFile, p.Object.AverageScans, p.Object.DeconvolutionTolerancePpm);
                }

                using (StreamWriter output = new StreamWriter(@"DeconvolutionOutput-" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture) + ".tsv"))
                {
                    output.WriteLine("Mass\tScore\tNumPeaks\tNumScans\tMinScan\tMaxScan\tElutionStart\tElutionEnd\tTotalNormalizedIntensity\tObservedCharges\tMostIntenseElutionTime\tMostIntenseCharge\tMostIntenseMz\tNumPeaksInMostIntenseEnvelope\tMostIntenseEnvelopeIntensity\tElutionOfMostIntenseCharge");

                    foreach (var nice in myMsDataFile.Deconvolute(p.Object.MinScan, p.Object.MaxScan, p.Object.MinAssumedChargeState, p.Object.MaxAssumedChargeState, p.Object.DeconvolutionTolerancePpm, p.Object.IntensityRatioLimit, p.Object.AggregationTolerancePpm, b => b.MsnOrder == 1).OrderByDescending(b => b.Score))
                    {
                        if ((nice.MaxScanIndex - nice.MinScanIndex + 1) >= p.Object.NumScansRequired)
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
        public int AverageScans { get; set; } = 1;
        public string FilePath { get; set; }
        public int NumScansRequired { get; set; } = 2;

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
            sb.AppendLine("AverageScans: " + AverageScans);
            sb.AppendLine("NumScansRequired: " + NumScansRequired);
            return sb.ToString();
        }
    }
}