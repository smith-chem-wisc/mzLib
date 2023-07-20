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
        public readonly ExtremumType ExtremumType;
        public Extremum(double retentionTime, int intensity, ExtremumType type)
        {
            RetentionTime = retentionTime;
            Intensity = intensity;
            ExtremumType = type;
        }

        public int CompareTo(Extremum other)
        {
            if (other == null) return 1;

            return RetentionTime.CompareTo(other.RetentionTime);
        }

        /// <summary>
        /// Returns true if both objects are extrema, both extrema 
        /// have the same intensity, same type, and retention times
        /// within 6 milliseconds of one another
        /// </summary>
        public override bool Equals(object obj)
        {
            return Equals(obj as Extremum);
        }
        //TODO: override GetHashcode

        /// <summary>
        /// Returns true if both extrema have the same intensity,
        /// same type, and retention times within 6 milliseconds
        /// of one another
        /// </summary>
        public bool Equals(Extremum obj)
        {
            return obj != null 
                && Math.Abs(obj - this) <= 0.0001 
                && obj.Intensity == this.Intensity
                && obj.ExtremumType == this.ExtremumType;
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
