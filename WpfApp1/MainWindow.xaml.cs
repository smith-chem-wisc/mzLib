using Chemistry;
using IO.ThermoRawFileReader;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using mzPlot;
using Nett;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Windows;

namespace WpfApp1
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            CreatePlot();
        }

        private void CreatePlot()
        {
            //var file = @"C:\Data\Yeast\02-15-17_YL-stnd_old-heat.raw";
            //var msFile = ThermoRawFileReader.LoadAllStaticData(file);


            //Plot p = new LinePlot(thePlotView, new List<Datum>());

            //double mass = 2965.35071;
            //for (int i = 0; i < 6; i++)
            //{
            //    double isotope = Constants.C13MinusC12 * i;
            //    var xic = msFile.ExtractIonChromatogram(mass + isotope, 3, new PpmTolerance(10), 132.18457);
            //    p.AddLinePlot(xic.Data, seriesTitle: (mass + isotope).ToString("F2"));
            //}

            Random r = new Random();

            //List<Datum> data = new List<Datum> {
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //    new Datum(r.Next(0, 1000), r.Next(0, 1000)),
            //};

            //var data = GetData().ToList();

            //Plot p = new LinePlot(thePlotView, data[0].Value, seriesTitle: data[0].Key);
            //p.AddScatterSeries(data[1].Value, seriesTitle: data[1].Key);
            //p.AddScatterSeries(data[2].Value, seriesTitle: data[2].Key);

            Stopwatch s = new Stopwatch();
            s.Start();

            StudentT n = new StudentT(0, 0.5, 10);
            List<Datum> data = new List<Datum>() 
            { 

            };
            for (double i = -1; i <= 1; i += 0.001)
            {
                for (double j = -1; j <= 1; j += 0.001)
                {
                    i = Math.Round(i, 3);
                    j = Math.Round(j, 3);

                    var density = n.Density(i) * n.Density(j) * 100;
                    data.Add(new Datum(i, j, density));
                }
            }

            //while (true)
            {
                //if (s.ElapsedMilliseconds > 1000)
                //{
                var p = new SurfacePlot(thePlotView, data, theImage);
                //    SurfacePlot.AngleX++;
                //    if (SurfacePlot.AngleX > 360)
                //    {
                //        SurfacePlot.AngleX = 0;
                //    }
                //    s.Restart();
                //}
            }



            //Toml.WriteFile(p, Path.Combine(@"C:\Data\Yeast\2020-09-29-20-17-44\Task1-SearchTask\ScatterPlot.toml"), Plot.tomlConfig);

        }

        private Dictionary<string, List<Datum>> GetData()
        {
            string path = @"C:\Data\Yeast\2020-09-29-20-17-44\Task1-SearchTask\AllPSMs.psmtsv";

            Dictionary<string, List<Datum>> dataPerFile = new Dictionary<string, List<Datum>>();

            using (StreamReader sr = new StreamReader(path))
            {
                int lineN = 0;
                while (sr.Peek() >= 0)
                {
                    lineN++;
                    var line = sr.ReadLine();

                    if (lineN == 1)
                    {
                        continue;
                    }

                    var split = line.Split(new char[] { '\t' });
                    var file = split[0];

                    if (!dataPerFile.TryGetValue(file, out var dataList))
                    {
                        dataPerFile.Add(file, new List<Datum>());
                    }

                    double qValue = double.Parse(split[48]);
                    if (qValue > 0.01)
                    {
                        break;
                    }

                    if (!double.TryParse(split[8], out double x))
                    {
                        continue;
                    }
                    if (!int.TryParse(split[11], out int notch) || notch != 0)
                    {
                        continue;
                    }
                    if (!double.TryParse(split[22], out double y))
                    {
                        continue;
                    }

                    dataPerFile[file].Add(new Datum(x, y));
                }
            }

            return dataPerFile;
        }
    }
}
