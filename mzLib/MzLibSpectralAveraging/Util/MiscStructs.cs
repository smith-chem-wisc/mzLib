using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibSpectralAveraging
{
    internal record struct BinnedPeak
    {
        internal BinnedPeak(int bin, double mz, double intensity, int specId)
        {
            Bin = bin;
            Mz = mz;
            Intensity = intensity;
            SpectraId = specId;
        }

        internal int Bin { get; init; }
        internal double Mz { get; init; }
        internal double Intensity { get; init; }
        internal int SpectraId { get; init; }
        public override string ToString()
        {
            return Mz + " : " + Intensity + " : " + SpectraId;
        }
    }
}
