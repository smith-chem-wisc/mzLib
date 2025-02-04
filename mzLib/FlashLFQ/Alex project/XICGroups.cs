using Easy.Common.Extensions;
using Microsoft.FSharp.Data.UnitSystems.SI.UnitNames;
using MzLibUtil;
using OpenMcdf;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using TopDownProteomics.IO.MzIdentMl;

namespace FlashLFQ.Alex_project
{
    public class XICGroups : IEnumerable<XICGroups>
    {
        public XIC reference;
        public List<XIC> XICs;
        public Dictionary<int, double> RTDict;
        public List<Extremum> sharedExtrema;
        public Dictionary<double, double> ExtremaInRef; // the shared extrema in the reference XIC (time/intensity)
        public List<Identification> XICGroup_IdList;
        public readonly Dictionary<double, Tuple<double,double>> indexedPeaks;

        /// <summary>
        /// Build the XIC groups with the reference XIC
        /// </summary>
        /// <param name="xICs"> The xICs list </param>
        /// <param name="extremaThreshold"> The minmun number to find the shared extrema </param>
        /// <param name="tolerance"> The tolerance window to pick up shared exterma </param>
        public XICGroups(List<XIC> xICs, double extremaThreshold = 0.55 , double tolerance = 0.10)
        {
            XICs = xICs;
            //this.ids = ids;
            reference = XICs.Where(p => p.Reference == true).First();    // set up the XIC reference
            RTDict = new Dictionary<int, double>();                      // build a dictionary to store the retention time shift of each XIC

            int xicID = 0;
            foreach (var xic in XICs)
            {
                RTDict.Add(xicID, xic.AlignXICs(reference)); //xic.AlignXICs(reference)
                xic.FindExtrema();
                xicID++;
            }

            sharedExtrema = new List<Extremum>();
            buildSharedExtrema(extremaThreshold, tolerance);    // find the sharedExtrema
            indexedPeaks = buildIndexedPeaks();                          // build the indexed peaks
            sharedExtremaInRef(reference, sharedExtrema);
            BuildIdList();
        }

        public void buildSharedExtrema(double count_threshold, double tolerance, double cutOff = 30000) 
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

            sharedExtrema = group_min.Concat(group_max).ToList();
            sharedExtrema.Sort((p1, p2) => p1.RetentionTime.CompareTo(p2.RetentionTime));
            Triming();
            int iii = 0;


        }

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

        private void Triming(double cutOff = 0.3) 
        {
            int removeIndex = 0;
            int count = sharedExtrema.Count();

            while (removeIndex < count - 1)
            {
                Extremum extre = sharedExtrema[removeIndex];
                Extremum extre_next = sharedExtrema[removeIndex + 1];

                double intensity = reference.SmoothedCubicSpline.Interpolate(extre.RetentionTime);
                double intensity_next = reference.SmoothedCubicSpline.Interpolate(extre_next.RetentionTime);
                double timeDiff = extre_next.RetentionTime - extre.RetentionTime;
                
                bool same_Min = extre.Type == ExtremumType.Minimum && extre.Type == extre_next.Type; // the constective minimum
                bool same_Max = extre.Type == ExtremumType.Maximum && extre.Type == extre_next.Type; // the constective minimum
                bool goUp = intensity_next > intensity; // the intensity is going up
                bool goDown = intensity_next < intensity; // the intensity is going down
                bool WithinTimeWindow = timeDiff < cutOff; // the time difference is within the cutOff


                if (same_Min && goUp && WithinTimeWindow)
                {
                    sharedExtrema.RemoveAt(removeIndex);
                    count--;
                }
                else if (same_Max && goDown && WithinTimeWindow) 
                {
                    sharedExtrema.RemoveAt(removeIndex);
                    count--;
                }
                else
                {
                    removeIndex++;
                }
            }

            int iii = 0;
            //while (removeIndex < count - 1)
            //{
            //    bool sameType = sharedExtrema[removeIndex].Type == sharedExtrema[removeIndex + 1].Type;
            //    double timeDiff = sharedExtrema[removeIndex + 1].RetentionTime - sharedExtrema[removeIndex].RetentionTime;
            //    if (sameType && timeDiff < cutOff)
            //    {
            //        sharedExtrema.RemoveAt(removeIndex + 1);
            //        count--;
            //    }
            //    else
            //    {
            //        removeIndex++;
            //    }
            //}
        }

        //public void buildSharedExtrema(double count_threshold , double tolerance , double cutOff = 0.16)
        //{
        //    List<Extremum> extremaList = new List<Extremum>();
        //    foreach (var xic in XICs)
        //    {
        //        extremaList.AddRange(xic.Extrema);
        //    }
        //    extremaList.Sort((p1,p2) => p1.RetentionTime.CompareTo(p2.RetentionTime));


        //    int index = 0;
        //    Dictionary<Extremum, List<Extremum>> group = new Dictionary<Extremum, List<Extremum>>();
        //    Extremum currentExtremum ;
        //    while (index < extremaList.Count()-1)
        //    {
        //        currentExtremum = extremaList[index]; 
        //        List<Extremum> currentGroup = new List<Extremum>() { currentExtremum };

        //        for (int i = index+1;  i < extremaList.Count() ; i++)
        //        {
        //            double timeDiff = extremaList[i].RetentionTime - currentExtremum.RetentionTime;
        //            index = i;

        //            if (timeDiff > tolerance)
        //            { 
        //                break;
        //            }

        //            else
        //            {
        //                if (extremaList[i].Type == currentExtremum.Type)
        //                {
        //                    currentGroup.Add(extremaList[i]);
        //                }
        //            }
        //        }

        //        if (group.ContainsKey(currentExtremum))
        //        {
        //            continue;
        //        }

        //        group.Add(currentExtremum, currentGroup);

        //    }

        //    sharedExtrema = group.Where(p => p.Value.Count() >= count_threshold * XICs.Count())
        //                         .Select(p=>p.Key).ToList();

        //    int removeIndex = 0;
        //    int count = sharedExtrema.Count();
        //    while (removeIndex < count-1) 
        //    {
        //        bool sameType = sharedExtrema[removeIndex].Type == sharedExtrema[removeIndex + 1].Type;
        //        double timeDiff = sharedExtrema[removeIndex + 1].RetentionTime - sharedExtrema[removeIndex].RetentionTime;
        //        if (sameType && timeDiff < cutOff)
        //        {
        //            sharedExtrema.RemoveAt(removeIndex + 1);
        //            count--;
        //        }
        //        else
        //        {
        //            removeIndex++;
        //        }
        //    }

        //    int I = 0;


        //}

        public void sharedExtremaInRef(XIC reference, List<Extremum> sharedExtre)
        {
            ExtremaInRef = new Dictionary<double, double>();
            foreach (var extre in sharedExtre)
            {
                double extreTime = extre.RetentionTime;
                double extreIntensity = reference.SmoothedCubicSpline.Interpolate(extreTime);
                ExtremaInRef.Add(extreTime, extreIntensity);
            }
        }

        /// <summary>
        /// Divide the shared extrema into individual peaks as indexed peaks
        /// </summary>
        /// <returns></returns>
        public Dictionary<double, Tuple<double, double>> buildIndexedPeaks(double windowLimit = 0.3)
        {
            var indexedPeaks = new Dictionary<double, Tuple<double, double>>();
            foreach (var extremPoint in sharedExtrema.Where(p=>p.Type == ExtremumType.Maximum))
            {
                double startTime = sharedExtrema.First().RetentionTime - 1;
                double endTime = sharedExtrema.Last().RetentionTime + 1;        
                int index_Peak = sharedExtrema.IndexOf(extremPoint);

                if (index_Peak - 1 >= 0 && sharedExtrema[index_Peak - 1].Type == ExtremumType.Minimum) // if the previous one is minimum, the start time is the minimum
                {
                    int indexLeft = index_Peak - 1;
                    startTime = sharedExtrema[indexLeft].RetentionTime;
                }

                else if(index_Peak - 1 >= 0)                                                          // if the previous one is maximum, the start time is the the midpint of the two maximum
                {
                    startTime = (sharedExtrema[index_Peak].RetentionTime + sharedExtrema[index_Peak - 1].RetentionTime)/2;
                }


                if (index_Peak + 1 < sharedExtrema.Count() && sharedExtrema[index_Peak + 1].Type == ExtremumType.Minimum) // if the next one is minimum, the end time is the minimum
                {
                    int indexRight = index_Peak + 1;
                    endTime = sharedExtrema[indexRight].RetentionTime;
                }

                else if(index_Peak + 1 < sharedExtrema.Count())                                                          // if the next one is maximum, the end time is the the midpint of the two maximum
                {
                    endTime = (sharedExtrema[index_Peak].RetentionTime + sharedExtrema[index_Peak + 1].RetentionTime)/2;
                }


                indexedPeaks.Add(extremPoint.RetentionTime, new Tuple<double, double> (startTime, endTime));
            }
            
            List<double> removeKeys = new List<double>();
            foreach (var window in indexedPeaks) 
            {
                double windowGap = window.Value.Item2 - window.Value.Item1;
                if (windowGap < windowLimit) 
                {
                    removeKeys.Add(window.Key);
                }
            }

            foreach (var key in removeKeys)
            {
                indexedPeaks.Remove(key);
            }

            return indexedPeaks;
        }


        public void BuildIdList() 
        {
            XICGroup_IdList = new List<Identification>();
            foreach (var xic in XICs) 
            {
                if (xic.Ids != null) 
                {
                    foreach (var id in xic.Ids) // Collect and modified the identifications in the XIC
                    {
                        MovedIdentification moveId = new MovedIdentification(id.FileInfo, id.BaseSequence, id.ModifiedSequence,
                            id.MonoisotopicMass, id.Ms2RetentionTimeInMinutes, id.PrecursorChargeState, id.ProteinGroups.ToList(), xic.RtShift);
                        moveId.PeakfindingMass = id.PeakfindingMass;
                        XICGroup_IdList.Add(moveId);
                    }
                }
            }
        }


        public IEnumerator<XICGroups> GetEnumerator()
        {
            throw new NotImplementedException();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            throw new NotImplementedException();
        }


    }
}
