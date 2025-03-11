using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ.Alex_project
{
    public class XICGroups : IEnumerable<XIC>
    {
        public XIC ReferenceXIC;
        public List<XIC> XICs;
        public Dictionary<int, double> RTDict;
        public List<Extremum> SharedExtrema;
        public Dictionary<double, double> ExtremaInRef; // the shared extrema in the reference XIC (time/intensity)
        public List<Identification> IdList;
        public readonly Dictionary<double, Tuple<double,double>> SharedPeaks; //format: retention time < start time, end time>

        /// <summary>
        /// Build the XIC groups with the reference XIC
        /// </summary>
        /// <param name="xICs"> The xICs list </param>
        /// <param name="sharedPeakThreshold"> The require number ratio to define the shared peaks, if threshold is 0.5, at lease 50% need to be tracked </param>
        /// <param name="tolerance"> The tolerance window to pick up shared exterma </param>
        public XICGroups(List<XIC> xICs, double sharedPeakThreshold = 0.55 , double tolerance = 0.10, double cutoff = 30000)
        {
            XICs = xICs;
            //this.ids = ids;
            ReferenceXIC = XICs.Where(p => p.Reference == true).First();    // set up the XIC reference
            RTDict = new Dictionary<int, double>();                      // build a dictionary to store the retention time shift of each XIC

            int xicID = 0;
            foreach (var xic in XICs)
            {
                RTDict.Add(xicID, xic.AlignXICs(ReferenceXIC)); //xic.AlignXICs(reference)
                xic.FindExtrema();
                xicID++;
            }

            SharedExtrema = new List<Extremum>();
            FindSharedExtrema(sharedPeakThreshold, tolerance, cutoff);    // find the sharedExtrema
            ProjectExtremaInRef(ReferenceXIC, SharedExtrema);         // project the sharedExtrema in the reference XIC
            SharedPeaks = BuildSharedPeaks();                   // build the indexed peaks
            
            BuildIdList();
        }

        /// <summary>
        /// Iterate the extrema from each XIC and generate the shared extrema from the XICs, the criteria are below
        /// (1) the intensity of the extrema should be larger the cutOff and
        /// (2) the time difference between each run should be within the tolerance
        /// (3) the number of the shared extrema should be larger than the threshold * number of XICs
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

            shared_min.Sort((p1, p2) => p1.RetentionTime.CompareTo(p2.RetentionTime));
            shared_max.Sort((p1, p2) => p1.RetentionTime.CompareTo(p2.RetentionTime));

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

            int iii = 0;
            
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
        public Dictionary<double, Tuple<double, double>> BuildSharedPeaks(double windowLimit = 0.3)
        {
            var sharedPeaks = new Dictionary<double, Tuple<double, double>>();
            foreach (var extremPoint in SharedExtrema.Where(p=>p.Type == ExtremumType.Maximum))
            {
                double startTime = SharedExtrema.First().RetentionTime - 1;
                double endTime = SharedExtrema.Last().RetentionTime + 1;        
                int maxPointIndex = SharedExtrema.IndexOf(extremPoint);

                if (maxPointIndex - 1 >= 0 && SharedExtrema[maxPointIndex - 1].Type == ExtremumType.Minimum) // if the previous one is minimum, the start time is the minimum
                {
                    int indexLeft = maxPointIndex - 1;
                    startTime = SharedExtrema[indexLeft].RetentionTime;
                }

                else if(maxPointIndex - 1 >= 0)                                                          // if the previous one is maximum, the start time is the the midpint of the two maximum
                {
                    startTime = (SharedExtrema[maxPointIndex].RetentionTime + SharedExtrema[maxPointIndex - 1].RetentionTime)/2;
                }


                if (maxPointIndex + 1 < SharedExtrema.Count() && SharedExtrema[maxPointIndex + 1].Type == ExtremumType.Minimum) // if the next one is minimum, the end time is the minimum
                {
                    int indexRight = maxPointIndex + 1;
                    endTime = SharedExtrema[indexRight].RetentionTime;
                }

                else if(maxPointIndex + 1 < SharedExtrema.Count())                                                          // if the next one is maximum, the end time is the the midpint of the two maximum
                {
                    endTime = (SharedExtrema[maxPointIndex].RetentionTime + SharedExtrema[maxPointIndex + 1].RetentionTime)/2;
                }


                sharedPeaks.Add(extremPoint.RetentionTime, new Tuple<double, double> (startTime, endTime));
            }

            // remove the peaks with the time window less than the window limit
            List<double> removedPeaks = new List<double>();
            foreach (var window in sharedPeaks) 
            {
                double windowGap = window.Value.Item2 - window.Value.Item1;
                if (windowGap < windowLimit) 
                {
                    removedPeaks.Add(window.Key);
                }
            }

            foreach (var removedPeak in removedPeaks)
            {
                sharedPeaks.Remove(removedPeak);
            }

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
                    foreach (var id in xic.Ids.Where(p=> p.FileInfo.Equals(xic.SpectraFile))) // Collect and modified the identifications in the XIC, Ignore the borrowed Id
                    {
                        // The moved id is the id with the modified retention time,
                        // ex. XIC no.1 (TimeShift is 0.5 min, id retentionTime is 10 mins), then the moved id retention time is 10.5 mins in reference XIC
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
