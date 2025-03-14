using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ.IsoTracker
{
    /// <summary>
    /// The XICGroups class is used to store the XICs the shared extrema, and the shared peak region between the XICs
    /// </summary>
    public class XICGroups : IEnumerable<XIC>
    {
        /// <summary>
        /// The reference XIC for alignment
        /// </summary>
        public XIC ReferenceXIC;
        /// <summary>
        /// The XICs list
        /// </summary>
        public List<XIC> XICs;
        /// <summary>
        /// The retention time shift of each XIC
        /// </summary>
        public Dictionary<int, double> RTDict;
        /// <summary>
        /// The shared extrema in the XICs, time project to the reference XIC
        /// </summary>
        public List<Extremum> SharedExtrema;
        /// <summary>
        /// The shared extrema in the reference XIC, we project the shared extrema in the reference XIC
        /// </summary>
        public Dictionary<double, double> ExtremaInRef;
        public List<Identification> IdList;
        /// <summary>
        /// The shared peak region in the XICGroups
        /// </summary>
        public List<PeakRegion> SharedPeaks { get; private set; }

        /// <summary>
        /// Build the XIC groups with the reference XIC
        /// </summary>
        /// <param name="xICs"> The xICs list </param>
        /// <param name="sharedPeakThreshold"> The require number ratio to define the shared peaks, if threshold is 0.5, at lease 50% need to be tracked </param>
        /// <param name="tolerance"> The tolerance window to pick up shared exterma </param>
        public XICGroups(List<XIC> xICs, double sharedPeakThreshold = 0.55 , double tolerance = 0.10, double cutOff = 30000)
        {
            XICs = xICs;
            //this.ids = ids;
            ReferenceXIC = XICs.First(p => p.Reference);    // set up the XIC reference
            RTDict = new Dictionary<int, double>();                             // build a dictionary to store the retention time shift of each XIC

            int xicID = 0;
            foreach (var xic in XICs)
            {
                RTDict.Add(xicID, xic.AlignXICs(ReferenceXIC)); //xic.AlignXICs(reference)
                xic.FindExtrema();
                xicID++;
            }

            SharedExtrema = new List<Extremum>();
            FindSharedExtrema(sharedPeakThreshold, tolerance, cutOff);    // find the sharedExtrema
            ProjectExtremaInRef(ReferenceXIC, SharedExtrema);         // project the sharedExtrema in the reference XIC
            SharedPeaks = BuildSharedPeaks();                   // build the indexed peaks
            
            BuildIdList();
        }

        /// <summary>
        /// Iterate the extrema from each XIC and generate the shared extrema from the XICs, the criteria are below
        /// (1) The intensity of the extrema should be larger the cutOff
        /// (2) The time difference between each run should be within the tolerance
        /// (3) The number of the shared extrema should be larger than the threshold * number of XICs
        /// </summary>
        /// <param name="count_threshold"> the minimum number for share tracking </param>
        /// <param name="tolerance"> Time difference tolerance</param>
        /// <param name="cutOff"> Intensity cutoff</param>
        public void FindSharedExtrema(double count_threshold, double tolerance, double cutOff) 
        {
            List<Extremum> shared_min = new List<Extremum>();
            List<Extremum> shared_max = new List<Extremum>();
            foreach (var xic in XICs)
            {
                shared_min.AddRange(xic.Extrema.Where(p => p.Type == ExtremumType.Minimum && p.Intensity >= cutOff));
                shared_max.AddRange(xic.Extrema.Where(p => p.Type == ExtremumType.Maximum && p.Intensity >= cutOff));
            }

            shared_min.Sort((p1, p2) => p1.CompareTo(p2));
            shared_max.Sort((p1, p2) => p1.CompareTo(p2));

            var group_min = SortExtrema(shared_min, tolerance, count_threshold);
            var group_max = SortExtrema(shared_max, tolerance, count_threshold);

            SharedExtrema = group_min.Concat(group_max).ToList();
            SharedExtrema.Sort((p1, p2) => p1.RetentionTime.CompareTo(p2.RetentionTime)); // sort the shared extrema by the retention time
            SharePeakTrimming();


        }

        /// <summary>
        /// Sort the extrema into groups based on the tolerance and the count threshold
        /// </summary>
        /// <param name="extremaList"></param>
        /// <param name="tolerance"></param>
        /// <param name="count_threshold"></param>
        /// <returns></returns>
        private List<Extremum> SortExtrema(List<Extremum> extremaList, double tolerance, double count_threshold) 
        {
            int index = 0;
            Dictionary<Extremum, List<Extremum>> group = new Dictionary<Extremum, List<Extremum>>();
            Extremum currentExtremum;
            while (index < extremaList.Count() - 1)
            {
                currentExtremum = extremaList[index];
                List<Extremum> currentGroup = new List<Extremum>() { currentExtremum };

                for (int i = index + 1; i < extremaList.Count(); i++)
                {
                    double timeDiff = extremaList[i].RetentionTime - currentExtremum.RetentionTime;
                    index = i;

                    if (timeDiff <= tolerance)
                    {
                        currentGroup.Add(extremaList[i]); ;
                    }

                    else
                    {
                        break;
                    }
                }

                if (group.ContainsKey(currentExtremum))
                {
                    continue;
                }

                group.Add(currentExtremum, currentGroup);

            }

            return group.Where(p => p.Value.Count() >= count_threshold * XICs.Count())
                                 .Select(p => p.Key).ToList();
        }

        /// <summary>
        /// Trim the shared extremas in the close time range. Ex. if the two minimums are close and the intensity is going up, remove the first one
        /// </summary>
        /// <param name="cutOff"> The time range for trimming</param>
        private void SharePeakTrimming(double cutOff = 0.3) 
        {
            int removeIndex = 0;
            int count = SharedExtrema.Count();

            while (removeIndex < count - 1)
            {
                Extremum extre = SharedExtrema[removeIndex];
                Extremum extre_next = SharedExtrema[removeIndex + 1];

                double intensity = ReferenceXIC.SmoothedCubicSpline.Interpolate(extre.RetentionTime);
                double intensity_next = ReferenceXIC.SmoothedCubicSpline.Interpolate(extre_next.RetentionTime);
                double timeDiff = extre_next.RetentionTime - extre.RetentionTime;
                
                bool same_Min = extre.Type == ExtremumType.Minimum && extre.Type == extre_next.Type; // the constective minimum
                bool same_Max = extre.Type == ExtremumType.Maximum && extre.Type == extre_next.Type; // the constective minimum
                bool goUp = intensity_next > intensity; // the intensity is going up
                bool goDown = intensity_next < intensity; // the intensity is going down
                bool WithinTimeWindow = timeDiff < cutOff; // the time difference is within the cutOff


                if (same_Min && goUp && WithinTimeWindow)
                {
                    SharedExtrema.RemoveAt(removeIndex);
                    count--;
                }
                else if (same_Max && goDown && WithinTimeWindow) 
                {
                    SharedExtrema.RemoveAt(removeIndex);
                    count--;
                }
                else
                {
                    removeIndex++;
                }
            }
        }

        /// <summary>
        /// Generate the shared extrema information (retention and intensity) in the reference XIC
        /// </summary>
        /// <param name="reference"></param>
        /// <param name="sharedExtre"></param>
        public void ProjectExtremaInRef(XIC reference, List<Extremum> sharedExtrema)
        {
            ExtremaInRef = new Dictionary<double, double>();
            foreach (var extre in sharedExtrema)
            {
                double extreTime = extre.RetentionTime;
                double extreIntensity = reference.SmoothedCubicSpline.Interpolate(extreTime);
                ExtremaInRef.Add(extreTime, extreIntensity);
            }
        }

        /// <summary>
        /// Use the shared extrema to build the shared peaks, group the peaks with consecutive minimum , maximum and minimum.
        /// </summary>
        /// /// <param name="windowLimit"> The minimum time window to generate a peak region</param>
        /// <returns></returns>
        public List<PeakRegion> BuildSharedPeaks(double windowLimit = 0.3)
        {
            var sharedPeaks = new List<PeakRegion>();
            foreach (var extremumPoint in SharedExtrema.Where(p=>p.Type == ExtremumType.Maximum))
            {
                double startTime = SharedExtrema.First().RetentionTime - 1;
                double endTime = SharedExtrema.Last().RetentionTime + 1;        
                int maxPointIndex = SharedExtrema.IndexOf(extremumPoint);

                // if the previous one is minimum, the start time is the minimum
                if (maxPointIndex - 1 >= 0 && SharedExtrema[maxPointIndex - 1].Type == ExtremumType.Minimum) 
                {
                    int indexLeft = maxPointIndex - 1;
                    startTime = SharedExtrema[indexLeft].RetentionTime;
                }

                // if the previous one is maximum, the start time is the midpoint of the two maximum
                else if (maxPointIndex - 1 >= 0)                                                          
                {
                    startTime = (SharedExtrema[maxPointIndex].RetentionTime + SharedExtrema[maxPointIndex - 1].RetentionTime)/2;
                }

                // if the next one is minimum, the end time is the minimum
                if (maxPointIndex + 1 < SharedExtrema.Count() && SharedExtrema[maxPointIndex + 1].Type == ExtremumType.Minimum) 
                {
                    int indexRight = maxPointIndex + 1;
                    endTime = SharedExtrema[indexRight].RetentionTime;
                }

                // if the next one is maximum, the end time is the midpoint of the two maximum
                else if (maxPointIndex + 1 < SharedExtrema.Count())                                                          
                {
                    endTime = (SharedExtrema[maxPointIndex].RetentionTime + SharedExtrema[maxPointIndex + 1].RetentionTime)/2;
                }

                sharedPeaks.Add(new PeakRegion(extremumPoint.RetentionTime, startTime, endTime));
            }

            // remove the peaks with the time window less than the window limit
            sharedPeaks.RemoveAll(p => p.Width < windowLimit);
            return sharedPeaks;
        }

        /// <summary>
        /// Because IsoTracker need the id from each runs for id sharing, we need to build the id list for each XIC group
        /// Note: in order to compare the id, we need to modify the retention time of the id based on the XIC shift. The timeline is Reference XIC
        /// </summary>
        public void BuildIdList() 
        {
            IdList = new List<Identification>();
            foreach (var xic in XICs) 
            {
                if (xic.Ids != null) 
                {
                    // Collect and modified the identifications in the XIC, Ignore the borrowed IDs
                    foreach (var id in xic.Ids.Where(p=> p.FileInfo.Equals(xic.SpectraFile))) 
                    {
                        // The moved id is the id with the modified retention time,
                        // ex. XIC no.1 (TimeShift is 0.5 min, id retentionTime is 10 min), then the moved id retention time is 10.5 min in reference XIC
                        MovedIdentification moveId = new MovedIdentification(id.FileInfo, id.BaseSequence, id.ModifiedSequence,
                            id.MonoisotopicMass, id.Ms2RetentionTimeInMinutes, id.PrecursorChargeState, id.ProteinGroups.ToList(), xic.RtShift);
                        moveId.PeakfindingMass = id.PeakfindingMass;
                        IdList.Add(moveId);
                    }
                }
            }
        }

        public IEnumerator<XIC> GetEnumerator() => XICs.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
    }
}
