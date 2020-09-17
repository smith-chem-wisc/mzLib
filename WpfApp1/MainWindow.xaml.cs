using BayesianEstimation;
using MathNet.Numerics.Distributions;
using mzPlot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MzLibUtil;
using MassSpectrometry;
using IO.MzML;
using IO.ThermoRawFileReader;
using Chemistry;

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
            ShowPlot();
        }

        public void ShowPlot()
        {
            //Datum dp1 = new Datum(1, 145, label: "test1");
            //Datum dp2 = new Datum(2, 200, label: "test2");
            //Datum dp3 = new Datum(5, 45, label: "test3");

            //Normal n = new Normal();
            //List<Datum> data = new List<Datum>();
            //for (int i = 0; i < 100; i++)
            //{
            //    data.Add(new Datum(n.Sample()));
            //}

            var filePath = ThermoRawFileReader.LoadAllStaticData(@"C:\Data\Yeast\02-15-17_YL-stnd_old-heat.raw");
            Plot s = new XicPlot(thePlotView, new ExtractedIonChromatogram(new List<Datum>()));

            double monoMass = 2965.35071;

            for (int i = 0; i < 5; i++)
            {
                var isotopeXic = Chromatography.ExtractIonChromatogram(
                    filePath, 
                    mass: monoMass + Constants.C13MinusC12 * i, 
                    charge: 3, 
                    massTolerance: new PpmTolerance(10), 
                    retentionTime: 132.18457);

                s.AddXicPlot(isotopeXic);
            }
        }
    }
}
