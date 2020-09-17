using MassSpectrometry;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class XicPlot : Plot
    {
        /// <summary>
        /// Creates an XIC plot (as a line plot).
        /// </summary>
        public XicPlot(PlotView oxyPlotView, ExtractedIonChromatogram xic) : base(oxyPlotView)
        {
            AddLinePlot(xic.Data);
        }
    }
}
