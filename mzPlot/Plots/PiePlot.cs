using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class PiePlot : Plot
    {
        /// <summary>
        /// Creates a pie plot. Not implemented yet.
        /// </summary>
        public PiePlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddPiePlot(data);
        }
    }
}
