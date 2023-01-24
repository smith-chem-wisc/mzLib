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

    public struct DoubleArray
    {
        public double[] Array { get; set; }

        public DoubleArray(double[] doubleArray)
        {
            Array = doubleArray;
        }

        public static double[] operator +(DoubleArray a, double[] b)
        {
            for (int i = 0; i < a.Array.Length; i++)
            {
                a.Array[i] += b[i]; 
            }
            return a.Array;
        }

        public static implicit operator DoubleArray(double[] array)
        {
            return new DoubleArray(array);
        }

        public static DoubleArray operator /(DoubleArray a, double scalar)
        {
            for (int i = 0; i < a.Array.Length; i++)
            {
                a.Array[i] /= scalar; 
            }

            return a; 
        }
        public static DoubleArray operator /(DoubleArray a, int scalar)
        {
            return a / (double)scalar; 
        }
    }
}
