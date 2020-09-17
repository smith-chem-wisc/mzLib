using MassSpectrometry;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class XicPlot : Plot
    {
        public XicPlot(PlotView oxyPlotView, ExtractedIonChromatogram xic) : base(oxyPlotView)
        {
            AddLinePlot(xic.Data);
        }
    }
}
