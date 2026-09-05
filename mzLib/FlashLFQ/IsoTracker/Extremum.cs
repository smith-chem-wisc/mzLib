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
            if (others == null)
            {
                return false;
            }

            if (this.Intensity == others.Intensity && Math.Abs(this - others) < 0.006 && this.Type == others.Type)
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Defers to the typed overload. The previous body was Equals(this, (Extremum)obj), which binds
        /// to the static object.Equals(object, object); that ends in a virtual objA.Equals(objB) and
        /// arrives straight back here, so any non-generic comparison of two distinct Extrema recursed
        /// until the stack ran out. Generic callers did not recurse, because Dictionary and LINQ resolve
        /// IEquatable&lt;Extremum&gt; and reach the typed overload directly -- but see GetHashCode: until it
        /// existed they were broken a different way, bucketing by reference so equal Extrema never met.
        /// </summary>
        public override bool Equals(Object obj)
        {
            return obj is Extremum other && Equals(other);
        }

        /// <summary>
        /// Hashes the two fields equality compares exactly. RetentionTime is excluded deliberately: it is
        /// compared within 0.006, and no hash can put two values that near each other in one bucket without
        /// putting everything there. Excluding it keeps the contract -- equal Extrema hash equally -- at the
        /// cost of collisions between same-intensity, same-type Extrema at unrelated retention times, which
        /// the typed Equals then separates.
        ///
        /// The tolerance still makes equality non-transitive (A and B within 0.006, B and C within 0.006, A
        /// and C not), so which Extrema a Dictionary treats as one key depends on insertion order. That is a
        /// property of the equality, not of the hash, and it is why Extremum is a poor dictionary key at all.
        /// </summary>
        public override int GetHashCode()
        {
            return HashCode.Combine(Intensity, Type);
        }

        public static double operator - (Extremum extremun1, Extremum extrenum2)
        {
            double rtDiff = extremun1.RetentionTime - extrenum2.RetentionTime;
            return rtDiff;
        }
    }
}
