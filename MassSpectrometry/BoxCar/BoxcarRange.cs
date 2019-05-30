using System;
using System.Collections;


namespace BoxCar
{
    /// <summary>
    /// A BoxcarRange has a Start and End. These are m/z values.
    /// This object is used for storing information about the boxcar scanning method.
    /// The information usually comes from the metadata about a scan.
    /// Each BoxcarRange can be stored in a SetOfBoxcarRanges object.
    /// </summary>
    public class BoxcarRange : IComparable
    {
        public double Start { get; set; }
        public double End { get; set; }

        /// <summary>
        /// Constructor for the Boxcar class.
        /// </summary>
        /// <param name="startpoint"></param>
        /// <param name="endpoint"></param>
        public BoxcarRange(double startpoint, double endpoint)
        {
            Start = startpoint;
            End = endpoint;
        }

        /// <summary>
        /// Default constructor for the Boxcar class.
        /// </summary>
        public BoxcarRange()
        {

        }

        // Methods/classes for sorting boxcars in ascending order of startpoint (NOT CURRENTLY USED)

        /// <summary>
        /// Compares the current boxcar object's start point to that of the input boxcar object.
        /// Returns 1 if the current boxcar object's start point is greater, -1 if the current boxcar's startpoint is less, and 0 if they are equal.
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        int IComparable.CompareTo(object obj)
        {
            BoxcarRange b = (BoxcarRange)obj;
            return String.Compare(this.Start.ToString(), b.Start.ToString());
        }

        /// <summary>
        /// helper class for the SortAscending() method
        /// </summary>
        private class SortAscendingHelper : IComparer
        {
            public int Compare(object x, object y)
            {
                throw new NotImplementedException();
            }

            int IComparer.Compare(object a, object b)
            {
                BoxcarRange b1 = (BoxcarRange)a;
                BoxcarRange b2 = (BoxcarRange)b;
                if (b1.Start > b2.Start)
                    return 1;
                if (b1.Start < b2.Start)
                    return -1;
                else
                    return 0;
            }
        }

        /// <summary>
        /// returns a helper object for the sort ascending method
        /// </summary>
        /// <returns></returns>
        public static IComparer SortAscending()
        {
            return (IComparer)new SortAscendingHelper();
        }
    }
}
