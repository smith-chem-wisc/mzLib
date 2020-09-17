using BayesianEstimation;
using MassSpectrometry;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class SpectrumPlot : Plot
    {
        public SpectrumPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddSpectrumPlot(data);
        }

        public SpectrumPlot(PlotView oxyPlotView, MzSpectrum spectrum, List<Datum> dataToAnnotate = null) : base(oxyPlotView)
        {
            List<Datum> spectrumData = new List<Datum>();

            for (int i = 0; i < spectrum.XArray.Length; i++)
            {
                double mz = spectrum.XArray[i];
                double intensity = spectrum.YArray[i];

                spectrumData.Add(new Datum(mz, intensity));
            }

            AddSpectrumPlot(spectrumData);

            // TODO: add annotation
        }
    }
}
