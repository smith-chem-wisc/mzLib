using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public abstract class MzPeak : Peak, IMzPeak
    {
        public MzPeak(double mz, double intensity)
            :base(mz, intensity)
        {
        }

        public double Intensity
        {
            get
            {
                return Y;
            }
        }

        public double Mz
        {
            get
            {
                return X;
            }
        }

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", X, Y);
        }

    }
}
