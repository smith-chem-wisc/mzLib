using System.Linq;
using MzLibUtil;

namespace SpectralAveraging
{
    public static class SpectraNormalization
    {
        /// <summary>
        ///     Calls specific normalization function
        /// </summary>
        /// <param name="yArrays">yArrays to be normalized</param>
        /// <param name="normalizationType">how to normalize spectra</param>
        /// <exception cref="MzLibException"></exception>
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

        /// <summary>
        ///     Divides each y value by its Tic value, sum of all y from one spectra will equal one
        /// </summary>
        /// <param name="yArrays">y arrays to be normalized</param>
        private static void NormalizeAbsoluteToTic(double[][] yArrays)
        {
            foreach (var t in yArrays)
            {
                var totalIonCurrent = t.Sum();
                if (totalIonCurrent == 0) totalIonCurrent = 1;

                for (var j = 0; j < t.Length; j++) t[j] = t[j] / totalIonCurrent;
            }
        }

        /// <summary>
        ///     Divides each y value by its Tic value, and multiples by the average Tic value
        /// </summary>
        /// <param name="yArrays">y arrays to be normalized</param>
        private static void NormalizeRelativeToTics(double[][] yArrays)
        {
            var tics = yArrays.Select(p => p.Sum()).ToArray();
            var averageTic = tics.Average();

            for (var i = 0; i < yArrays.Length; i++)
            {
                var tic = tics[i] == 0 ? 1 : tics[i];
                for (var j = 0; j < yArrays[i].Length; j++) yArrays[i][j] = yArrays[i][j] / tic * averageTic;
            }
        }
    }
}