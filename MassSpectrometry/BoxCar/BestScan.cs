using MassSpectrometry;

/// <summary>
/// Written by Nicole Frey, May-June 2019 for the Smith Group in the UW Madison chemistry department, with direction from Leah Schaffer.
/// </summary>
namespace BoxCar
{
    /// <summary>
    /// A BestScan object holds information about a scan that matches up with a BoxcarRange object.
    /// 
    /// How this is used in the program:
    /// The FirstIndex of a BestScan is the position of the first m/z value (in the x array of that scan) which is greater than the startpoint of a BoxcarRange object.
    /// The LastIndex of a BestScan is the position of the last m/z value (in the x array of that scan) which is less than or equal to the startpoint of a BoxcarRange object.
    /// TotalIntensity is the combined intensity (from the y array of that scan) corresponding to all the m/z values from the FirstIndex to the LastIndex (inclusive).
    /// </summary>
    public class BestScan
    {
        public int FirstIndex { get; set; }
        public int LastIndex { get; set; }
        public double TotalIntensity { get; set; }
        public MsDataScan Scan { get; set; }

        /// <summary>
        /// Constructor for the BestScan class
        /// </summary>
        /// <param name="firstIndex"></param>
        /// <param name="lastIndex"></param>
        /// <param name="totalIntensity"></param>
        /// <param name="scan"></param>
        public BestScan(int firstIndex, int lastIndex, double totalIntensity, MsDataScan scan)
        {
            FirstIndex = firstIndex;
            LastIndex = lastIndex;
            TotalIntensity = totalIntensity;
            Scan = scan;
        }
    }
}
