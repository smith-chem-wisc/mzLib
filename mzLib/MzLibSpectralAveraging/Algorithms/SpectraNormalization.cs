using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MzLibSpectralAveraging
{
    public static class SpectraNormalization
    {
        public static void NormalizeSpectra(double[][] yArrays, NormalizationType normalizationType)
        {
            switch (normalizationType)
            {
                case NormalizationType.NoNormalization:
                    return;

                case NormalizationType.AbsoluteToTic:
                    NormalizeAbsoluteToTic(yArrays);
                    break;

                case NormalizationType.RelativeToTics:
                    NormalizeRelativeToTics(yArrays);
                    break;

                default:
                    throw new MzLibException("Normalization Type not yet implemented");
            }
        }

        private static void NormalizeAbsoluteToTic(double[][] yArrays)
        {
            for (int i = 0; i < yArrays.Length; i++)
            {
                double totalIonCurrent = yArrays[i].Sum();
                if (totalIonCurrent == 0) totalIonCurrent = 1;

                for (int j = 0; j < yArrays[i].Length; j++)
                {
                    yArrays[i][j] = yArrays[i][j] / totalIonCurrent;
                }
            }
        }

        private static void NormalizeRelativeToTics(double[][] yArrays)
        {
            var tics = yArrays.Select(p => p.Sum()).ToArray();
            var averageTic = tics.Average();

            for (int i = 0; i < yArrays.Length; i++)
            {
                var tic = tics[i] == 0 ? 1 : tics[i];
                for (int j = 0; j < yArrays[i].Length; j++)
                {
                    yArrays[i][j] = yArrays[i][j] / tic * averageTic;
                }
            }
        }

    }
}
