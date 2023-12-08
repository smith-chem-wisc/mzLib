using MzLibUtil;
using System;

namespace Readers.Bruker.TimsTofReader
{
    public class TofPeak : IComparable
    {
        public double Mz { get; }
        public int Intensity { get; }

        public TofPeak(double mz, int intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public static TofPeak operator +(TofPeak peak1, TofPeak peak2)
            => new TofPeak((peak1.Mz + peak2.Mz) / 2, peak1.Intensity + peak2.Intensity);

        public int CompareTo(Object obj)
        {
            if (obj == null) return 1;

            TofPeak otherPeak = obj as TofPeak;
            if (otherPeak == null) throw new ArgumentException("Object is not a TofPeak");
            return this.Mz.CompareTo(otherPeak.Mz);
        }

        public int CompareTo(TofPeak otherPeak, Tolerance tol)
        {
            if (tol.Within(this.Mz, otherPeak.Mz)) return 0; // peaks are equal within tolerance
            return this.Mz.CompareTo(otherPeak.Mz);
        }

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", Mz, Intensity);
        }
    }
}