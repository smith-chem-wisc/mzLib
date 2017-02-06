using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Spectra
{
    public abstract class Peak : IPeak
    {
        public double X { get; protected set; }
        public double Y { get; protected set; }
    }
}
