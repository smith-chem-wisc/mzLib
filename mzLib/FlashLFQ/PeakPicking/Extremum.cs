using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace FlashLFQ.PeakPicking
{
    public class Extremum : IComparable<Extremum>
    {
        /// <summary>
        /// Aligned retention time
        /// </summary>
        public readonly double RetentionTime;
        /// <summary>
        /// Interpolated intensity
        /// </summary>
        public readonly int Intensity;
        public readonly ExtremumType Type;
        public Extremum(double retentionTime, int intensity, ExtremumType type)
        {
            RetentionTime = retentionTime;
            Intensity = intensity;
            Type = type;
        }

        public int CompareTo(Extremum other)
        {
            if (other == null) return 1;

            return RetentionTime.CompareTo(other.RetentionTime);
        }

        

        public static bool operator > (Extremum operand1, Extremum operand2)
        {
            return operand1.CompareTo(operand2) > 0;
        }

        public static bool operator < (Extremum operand1, Extremum operand2)
        {
            return operand1.CompareTo(operand2) < 0;
        }

        public static bool operator >= (Extremum operand1, Extremum operand2)
        {
            return operand1.CompareTo(operand2) >= 0;
        }

        public static bool operator <= (Extremum operand1, Extremum operand2)
        {
            return operand1.CompareTo(operand2) <= 0;
        }

        public static double operator - (Extremum operand1, Extremum operand2)
        {
            return operand1.RetentionTime - operand2.RetentionTime;
        }


    }

    public enum ExtremumType
    {
        Maximum,
        Minimum
    }
}
