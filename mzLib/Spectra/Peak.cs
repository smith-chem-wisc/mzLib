using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Spectra
{
    public abstract class Peak : IPeak
    {
        public double X { get; private set; }
        public double Y { get; private set; }
        public Peak(double x, double y)
        {
            X = x;
            Y = y;
        }
    }
}
