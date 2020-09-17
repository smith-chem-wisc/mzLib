using MzLibUtil;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class BoxPlot : Plot
    {
        public BoxPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddHistogram(data);
        }
    }
}
