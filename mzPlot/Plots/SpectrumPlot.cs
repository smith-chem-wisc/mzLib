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
        /// <summary>
        /// Plots a mass spectrum. The data X value is the m/z of the spectral line, the data Y value is the intensity (height) of the spectral line.
        /// </summary>
        public SpectrumPlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddSpectrumPlot(data);
        }

        public SpectrumPlot(PlotView oxyPlotView, MzSpectrum spectrum) : base(oxyPlotView)
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
