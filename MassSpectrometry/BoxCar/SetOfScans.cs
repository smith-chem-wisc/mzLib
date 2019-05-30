using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;


namespace BoxCar
{
    /// <summary>
    /// Each SetOfScans object has three lists containing MsDataScan objects.
    /// Ms1scans stores the full ms (ms1) scans in the set.
    /// BoxcarScans stores the scans in the set that were taken with the boxcar method.
    /// Ms2scans stores the ms2 scans in the set.
    /// 
    /// When the scan is taken, each SetOfScans object corresponds to a box, which can have multiple boxcars in it. 
    /// If the scan is preset to "12 boxes" or some other number of boxes, there will be an equal number of SetOfScans objects in the Program
    /// (unless some are incomplete, in which case they will be removed).
    /// </summary>
    public class SetOfScans : IEnumerable
    {
        public List<MsDataScan> Ms1scans { get; set; } // this will likely have 1-2 scans in it    
        public List<MsDataScan> BoxcarScans { get; set; } // 2-3 scans probably
        public List<MsDataScan> Ms2scans { get; set; } // many scans

        /// <summary>
        /// Constructor for the SetOfScans object with full ms1, boxcar, and ms2 scans.
        /// </summary>
        /// <param name="ms1Scans"></param>
        /// <param name="boxcarScans"></param>
        /// <param name="ms2Scans"></param>
        public SetOfScans(List<MsDataScan> ms1_list, List<MsDataScan> boxcar_list, List<MsDataScan> ms2_list)
        {
            Ms1scans = ms1_list;
            Ms2scans = ms2_list;
            BoxcarScans = boxcar_list;
        }

        /// <summary>
        /// Constructor for the SetOfScans object with full ms1 and boxcar scans only
        /// </summary>
        /// <param name="ms1_list"></param>
        /// <param name="boxcar_list"></param>
        public SetOfScans(List<MsDataScan> ms1_list, List<MsDataScan> boxcar_list)
        {
            Ms1scans = ms1_list;
            BoxcarScans = boxcar_list;
            Ms2scans = new List<MsDataScan>();
        }

        /// <summary>
        /// Default Constructor for the SetOfScans object
        /// </summary>
        public SetOfScans()
        {
            Ms1scans = new List<MsDataScan>();
            BoxcarScans = new List<MsDataScan>();
            Ms2scans = new List<MsDataScan>();
        }

        /// <summary>
        /// Adds an MsDataScan to the end of the Ms1scans list.
        /// </summary>
        /// <param name="scan"></param>
        public void AddToMs1Scans(MsDataScan scan)
        {
            Ms1scans.Add(scan);
        }

        /// <summary>
        /// Adds an MsDataScan to the end of the BoxcarScans list.
        /// </summary>
        /// <param name="scan"></param>
        public void AddToBoxcarScans(MsDataScan scan)
        {
            BoxcarScans.Add(scan);
        }
        
        /// <summary>
        /// Adds an MsDataScan to the end of the Ms2scans list.
        /// </summary>
        /// <param name="scan"></param>
        public void AddToMs2Scans(MsDataScan scan)
        {
            Ms2scans.Add(scan);
        }
        

        // Methods for iterating through the boxcars in a SetOfScans

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (IEnumerator)GetEnumerator();
        }

        public SetOfScansEnum GetEnumerator()
        {
            return new SetOfScansEnum(BoxcarScans);
        }

        public class SetOfScansEnum : IEnumerator
        {
            public List<MsDataScan> Boxcars;
            int position = -1;

            public SetOfScansEnum(List<MsDataScan> list)
            {
                Boxcars = list;
            }

            public bool MoveNext()
            {
                position++;
                return (position < Boxcars.Count);
            }

            public void Reset()
            {
                position = -1;
            }

            object IEnumerator.Current
            {
                get
                {
                    return Current;
                }
            }

            public MsDataScan Current
            {
                get
                {
                    try
                    {
                        return Boxcars.ElementAt(position);
                    }
                    catch (IndexOutOfRangeException)
                    {
                        throw new InvalidOperationException();
                    }
                }
            }
        }

        
    }
}
