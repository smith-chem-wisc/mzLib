using System;


namespace FlashLFQ.IsoTracker
{
    /// <summary>
    /// The type of the Extremum point, either a minimum or maximum
    /// </summary>
    public enum ExtremumType { Minimum, Maximum };

    /// <summary>
    /// An Extremum object stores the location of a local maxima or minima in an XIC
    /// </summary>
    public class Extremum : IComparable<Extremum>, IEquatable<Extremum>
    {
        /// <summary>
        /// The interpolated intensity of the Extremum point
        /// </summary>
        public readonly double Intensity;
        /// <summary>
        /// The retention time of the Extremum point
        /// </summary>
        public readonly double RetentionTime;
        /// <summary>
        /// The type of the Extremum point
        /// </summary>
        public readonly ExtremumType Type;

        public Extremum(double intensity, double retentionTime, ExtremumType type)
        {
            Intensity = intensity;
            RetentionTime = retentionTime;
            Type = type;
        }

        public int CompareTo(Extremum others)
        {
            if (others == null)
            {
                return 1;
            }

            if (this - others > 0)
            {
                return 1;
            }

            else if (this - others < 0)
            {
                return -1;
            }

            else 
            {
                return 0;
            }

        }

        public bool Equals(Extremum others)
        {
            if (this.Intensity == others.Intensity && Math.Abs(this - others) < 0.006 && this.Type == others.Type)
            {
                return true;
            }

            return false;
        }

        public override bool Equals(Object obj)
        {
            return Equals(this, (Extremum)obj);
        }

        public static double operator - (Extremum extremun1, Extremum extrenum2)
        {
            double rtDiff = extremun1.RetentionTime - extrenum2.RetentionTime;
            return rtDiff;
        }
    }
}
