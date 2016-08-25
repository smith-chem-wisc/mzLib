using Chemistry;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;

namespace Benchmark
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            string gitStatus = string.Empty;
            Stream stream = null;
            try
            {
                stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Benchmark." + "status.txt");
                using (StreamReader reader = new StreamReader(stream))
                {
                    stream = null;
                    // Use the reader object...
                    gitStatus = reader.ReadToEnd();
                }
            }
            finally
            {
                if (stream != null)
                    stream.Dispose();
            }
            string gitShow = string.Empty;
            stream = null;
            try
            {
                stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Benchmark." + "show.txt");
                using (StreamReader reader = new StreamReader(stream))
                {
                    stream = null;
                    // Use the reader object...
                    gitShow = reader.ReadToEnd();
                }
            }
            finally
            {
                if (stream != null)
                    stream.Dispose();
            }
            string compileTime = string.Empty;
            stream = null;
            try
            {
                stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Benchmark." + "buildDate.txt");
                using (StreamReader reader = new StreamReader(stream))
                {
                    stream = null;
                    // Use the reader object...
                    compileTime = reader.ReadToEnd();
                }
            }
            finally
            {
                if (stream != null)
                    stream.Dispose();
            }
            PopulatePeriodicTable();

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"..\..\..\Benchmark.txt"))
            {
                file.WriteLine("At compile time, date and time:");
                file.WriteLine(compileTime);
                file.WriteLine("");
                file.WriteLine("At compile time, git show was:");
                file.WriteLine(gitShow);
                file.WriteLine("");
                file.WriteLine("At compile time, git status was:");
                file.WriteLine(gitStatus);
                file.WriteLine("");

                BenchmarkFormula(file);
                file.WriteLine("");
                BenchmarkFormula2(file);
                file.WriteLine("");
                BenchmarkTimeGettingElementFromPeriodicTable(file);
                file.WriteLine("");
                BenchmarkGettingIsotopes(file);
                file.WriteLine("");
                BenchmarkIsotopicDistribution(file);
            }
        }

        private static void BenchmarkIsotopicDistribution(StreamWriter file)
        {
            file.WriteLine("Starting benchmark BenchmarkIsotopicDistribution");

            int numRepetitions = 100;

            Stopwatch stopWatch = new Stopwatch();

            var a = new ChemicalFormula("H100C100N100O100S100");
            double b = 0;
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                b += new IsotopicDistribution(a).Intensities.First();
            }
            stopWatch.Stop();
            file.WriteLine("Time for generating isotopic distributions: " + stopWatch.Elapsed + " a = " + a);

            file.WriteLine("Benchmark BenchmarkIsotopicDistribution finished");
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

        private static void BenchmarkGettingIsotopes(StreamWriter file)
        {
            file.WriteLine("Starting benchmark BenchmarkGettingIsotopes");

            int numRepetitions = 10000000;

            Stopwatch stopWatch = new Stopwatch();

            long a = 0;
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                a += PeriodicTable.GetElement(20).Isotopes.Count();
            }
            stopWatch.Stop();
            file.WriteLine("Time for getting isotopes1: " + stopWatch.Elapsed + " a = " + a);

            file.WriteLine("Benchmark BenchmarkGettingIsotopes finished");
        }

        private static void BenchmarkFormula(StreamWriter file)
        {
            file.WriteLine("Starting benchmark BenchmarkFormula");

            int numRepetitions = 100000;

            Stopwatch stopWatch = new Stopwatch();

            var a = new ChemicalFormula("H1H{1}10 H{2}10 O20 O{16}20 O{17}20 O{18}20 C{12}100 C100 C{13}100 S{32}200 S200 S{33}200 S{34}200 S{36}200");
            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var b = a.Formula + i;
            }
            stopWatch.Stop();
            file.WriteLine("Time for getting formulas: " + stopWatch.Elapsed);

            file.WriteLine("Benchmark BenchmarkFormula finished");
        }

        private static void BenchmarkFormula2(StreamWriter file)
        {
            file.WriteLine("Starting benchmark BenchmarkFormula2");

            int numRepetitions = 100000;

            Stopwatch stopWatch = new Stopwatch();

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = new ChemicalFormula("H" + i + "H{1}10 H{2}10 O20 O{16}20 O{17}20 O{18}20 C{12}100 C100 C{13}100 S{32}200 S200 S{33}200 S{34}200 S{36}200");
                var b = a.Formula + i;
            }
            stopWatch.Stop();
            file.WriteLine("Time for creating and getting formulas: " + stopWatch.Elapsed);

            file.WriteLine("Benchmark BenchmarkFormula2 finished");
        }

        private static void BenchmarkTimeGettingElementFromPeriodicTable(StreamWriter file)
        {
            file.WriteLine("Starting benchmark BenchmarkTimeGettingElementFromPeriodicTable");

            int numRepetitions = 100000000;

            Stopwatch stopWatch = new Stopwatch();

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = PeriodicTable.GetElement(1);
                var b = a.Protons + a.AverageMass + 4;
            }
            stopWatch.Stop();
            file.WriteLine("Time for getting by atomic number: " + stopWatch.Elapsed);

            stopWatch.Restart();
            for (int i = 0; i < numRepetitions; i++)
            {
                var a = PeriodicTable.GetElement("H");
                var b = a.Protons + a.AverageMass + 4;
            }
            stopWatch.Stop();
            file.WriteLine("Time for getting by atomic symbol: " + stopWatch.Elapsed);

            file.WriteLine("Benchmark BenchmarkTimeGettingElementFromPeriodicTable finished");
        }
    }
}