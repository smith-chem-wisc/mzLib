using Easy.Common.Extensions;
using SharpLearning.InputOutput.Csv;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Plotly.NET;
using Plotly.NET.CSharp;
using Chart = Plotly.NET.CSharp.Chart;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using System;
using System.Collections.Generic;
using System.Linq;
using Plotly.NET.LayoutObjects;
using Plotly.NET.TraceObjects;
using System.ComponentModel;
using System.Drawing;
using Color = Plotly.NET.Color;



namespace FlashLFQ.Alex_project
{
    public class ParallelSearch
    {

        readonly int threadNumber = Environment.ProcessorCount-1;
        readonly int[] threads;
        Dictionary<int, XIC> xicDict;
        XICGroups[] xicGroups;

        public ParallelSearch(Dictionary<int, XIC> xicDict)
        {
            threads = Enumerable.Range(0, threadNumber).ToArray();
            this.xicDict = xicDict;
            xicGroups = new XICGroups[threadNumber];
        }

        public void run()
        {
            Parallel.ForEach(threads, (currentThread) =>
            {
                List<XIC> xics = GroupedXIC(xicDict, currentThread);
                if (xics.Count() != 0)
                {
                    xicGroups[currentThread] = new XICGroups(xics, 0.5, 0.1);
                    draw(currentThread, xicGroups[currentThread]);
                }              
                
            });

            Console.WriteLine("The total sum is ");

        }

        /// <summary>
        /// Try to group the XICs from a Big XIC data set. If we met the reference XIC, we will build the XICGroups.
        /// </summary>
        /// <param name="xicDict"></param>
        public static List<XIC> GroupedXIC(Dictionary<int, XIC> xicDict, int thread)
        {
            Dictionary<int, XIC> xicsToGroup = new Dictionary<int, XIC>();

            lock (xicDict) 
            {
                foreach (var xic in xicDict)
                {
                    if (xic.Value.Reference == true && xicsToGroup.Where(p => p.Value.Reference).Count() > 0)
                    {
                        break;
                    }

                    else
                    {
                        xicsToGroup.Add(xic.Key, xic.Value);
                        xicDict.Remove(xic.Key);
                    }
                }
            }

            return xicsToGroup.Select(p => p.Value).ToList();

        }

        public static void draw(int currentThread, XICGroups xicGroup)
        {
            XIC referenceXIC = xicGroup.reference;

            var plotStack = Chart.Combine(
                xicGroup.XICs.Select(xic =>
                    Chart.Scatter<double, double, string>(
                        xic.Ms1Peaks.Select(p => p.RetentionTime).ToArray(),
                        xic.Ms1Peaks.Select(p => p.Intensity).ToArray(),
                        StyleParam.Mode.Lines_Markers,
                        Name: "No" + currentThread
                    )
                ).ToArray()

          )
          .WithTitle("XIC_differentTime")
          .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
          .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
          .WithSize(1200, 400);


            var plot_extrema = Chart.Combine(new[]
             {
                Chart.Scatter<double, double, string>(referenceXIC.Ms1Peaks.Select(P=>P.RetentionTime), referenceXIC.Ms1Peaks.Select(P=>P.Intensity),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, Name: "refernece"),
                Chart.Scatter<double, double, string>(xicGroup.ExtremaInRef.Select(p=> 
                p.Key), xicGroup.ExtremaInRef.Select(p=> p.Value),
                StyleParam.Mode.Markers, MarkerSymbol: StyleParam.MarkerSymbol.Star, Name: "shared Extrema").WithMarkerStyle(Size: 15),
            })

                 .WithTitle("XIC_sharedExtrema")
                 .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                 .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
                 .WithSize(1200, 400);


            var stack = Chart.Grid(new[] { plotStack, plot_extrema }, 2, 1)
                .WithSize(1200, 800);

            GenericChartExtensions.Show(stack);

        }

    }
}
