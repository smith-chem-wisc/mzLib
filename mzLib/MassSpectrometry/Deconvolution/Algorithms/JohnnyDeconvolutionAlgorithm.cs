using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    public class JohnnyDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
        : DeconvolutionAlgorithm(deconParameters)
    {
        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var firstIndex = spectrum.GetClosestPeakIndex(range.Minimum);
            var lastIndex = spectrum.GetClosestPeakIndex(range.Maximum);

            var mzs = spectrum.XArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();
            var intensities = spectrum.YArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();

            foreach (var johnnyDeconObjectType in ReplaceMeWithDllCall())
            {
                // pull johnny decon object type info out

                IsotopicEnvelope isotopicEnvelope = null;
                yield return isotopicEnvelope;
            }
        }

        private IEnumerable<object> ReplaceMeWithDllCall()
        {
            for (int i = 0; i < 10; i++)
            {
                yield return new object();
            }
        }
    }
}
