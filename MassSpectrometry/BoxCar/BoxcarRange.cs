using System;
using System.Collections;

/// <summary>
/// Written by Nicole Frey, May-June 2019 for the Smith Group in the UW Madison chemistry department, with direction from Leah Schaffer.
/// </summary>
namespace BoxCar
{
    /// <summary>
    /// A BoxcarRange has a Start and End. These are m/z values.
    /// This object is used for storing information about the boxcar scanning method.
    /// The information usually comes from the metadata about a scan.
    /// Each BoxcarRange can be stored in a SetOfBoxcarRanges object.
    /// </summary>
    public class BoxcarRange 
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
    }
}
