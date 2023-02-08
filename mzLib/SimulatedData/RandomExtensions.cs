using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.FSharp.Core;

namespace SimulatedData
{
    public static class RandomExtensions
    {
        public static double NextDouble(this Random random, double low, double high)
        {
            return random.NextDouble() * (high - low) + low;
        }
    }

    
}
