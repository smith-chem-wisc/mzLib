using Easy.Common.Extensions;
using MzLibUtil;
using OpenMcdf;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace FlashLFQ.Alex_project
{
    public class XICGroups : IEnumerable<XICGroups>
    {
        public XIC reference;
        public List<XIC> XICs;
        public Dictionary<int, double> RTDict;
        public List<Extremum> sharedExtrema;
        public Dictionary<double, double> ExtremaInRef; // the shared extrema in the reference XIC (time/intensity)
        readonly List<Identification> ids;
        public readonly Dictionary<double, Tuple<double,double>> indexedPeaks;

        /// <summary>
        /// Build the XIC groups with the reference XIC
        /// </summary>
        /// <param name="xICs"> The xICs list </param>
        /// <param name="extremaThreshold"> The minmun number to find the shared extrema </param>
        /// <param name="tolerance"> The tolerance window to pick up shared exterma </param>
        public XICGroups(List<XIC> xICs, double extremaThreshold = 0.5, double tolerance = 0.05)
        {
            XICs = xICs;
            //this.ids = ids;
            reference = XICs.Where(p => p.Reference == true).First();    // set up the XIC reference
            RTDict = new Dictionary<int, double>();                      // build a dictionary to store the retention time shift of each XIC

            int xicID = 0;
            foreach (var xic in XICs)
            {
                RTDict.Add(xicID, xic.AlignXICs(reference));
                xic.BuildSmoothedCubicSpline();
                xic.FindExtrema();
                xicID++;
            }

            sharedExtrema = new List<Extremum>();
            buildSharedExtrema_2(extremaThreshold, tolerance);    // find the sharedExtrema
            indexedPeaks = buildIndexedPeaks();                          // build the indexed peaks
            sharedExtremaInRef(reference, sharedExtrema);
        }

        /// <summary>
        /// Build the XICGroups, default extremaCutoff is 50% , default tolerance is 0.05 min
        /// </summary>
        /// <param name="xICs"></param>
        public XICGroups(List<XIC> xICs)
        {
            XICs = xICs;
            reference = XICs.Where(p => p.Reference == true).First();  
            RTDict = new Dictionary<int, double>();

            int xicID = 0;
            foreach (var xic in XICs)
            {
                RTDict.Add(xicID, xic.AlignXICs(reference));
                xic.BuildSmoothedCubicSpline();
                xic.FindExtrema();
                xicID++;
            }

            sharedExtrema = new List<Extremum>();
            buildSharedExtrema_2(0.5, 0.05); //default value, at least 50% of the XICs share the extrema
            indexedPeaks = buildIndexedPeaks();                          // build the indexed peaks
            sharedExtremaInRef(reference, sharedExtrema);
        }



        public void buildSharedExtrema(double count_threshold, double tolerance)
        {
           
            foreach (var ref_extremum in reference.Extrema)
            {
                
                int extremum_account = 0;

                foreach (var xic in XICs)
                {
                    foreach (var xic_extremum in xic.Extrema)
                    {
                        if (within(xic_extremum, ref_extremum, tolerance) && xic_extremum.Type == ref_extremum.Type)
                        {
                            extremum_account++;
                            break;
                        } 
                    }

                }

                if (extremum_account >= count_threshold * XICs.Count())  // if the extremum is shared by more than 50% of the XICs, then we accept
                {
                    sharedExtrema.Add(ref_extremum);
                }

            }
            
        }

        /// <summary>
        /// Check the extremum is within the local maxmun window
        /// </summary>
        /// <param name="item"> The extrema that want to compare </param>
        /// <param name="local_Max"> The local maximum in the reference </param>
        /// <returns> True: the extremun is within the window, False: the extremun is not within the window </returns>
        private bool within(Extremum item, Extremum local_Max, double tolerance = 0.05)
        {;

            double leftWindow = (local_Max.RetentionTime - tolerance >= 0)  ? local_Max.RetentionTime - tolerance : 0;


            double rightWindow = (local_Max.RetentionTime + tolerance <= reference.Ms1Peaks.Last().RetentionTime)
            ? local_Max.RetentionTime + tolerance : reference.Ms1Peaks.Last().RetentionTime;


            if (item.RetentionTime >= leftWindow && item.RetentionTime <= rightWindow && item.Type == local_Max.Type)
            {
                return true;
            }

            return false;
        }


        public void buildSharedExtrema_2(double count_threshold = 0.5, double tolerance = 0.1)
        {
            List<Extremum> extremaList = new List<Extremum>();
            foreach (var xic in XICs)
            {
                extremaList.AddRange(xic.Extrema);
            }
            extremaList.Sort((p1,p2) => p1.RetentionTime.CompareTo(p2.RetentionTime));


            int index = 0;
            Dictionary<Extremum, List<Extremum>> group = new Dictionary<Extremum, List<Extremum>>();

            while (index < extremaList.Count()-1)
            {
                var currentExtremum = extremaList[index]; 
                List<Extremum> currentGroup = new List<Extremum>();
                currentGroup.Add(currentExtremum);

                for (int i = index+1;  i < extremaList.Count() ; i++)
                {
                    double timeDiff = extremaList[i].RetentionTime - currentExtremum.RetentionTime;
                    index = i;

                    if (timeDiff > tolerance)
                    { 
                        break;
                    }

                    else
                    {
                        if (extremaList[i].Type == currentExtremum.Type)
                        {
                            currentGroup.Add(extremaList[i]);
                        }
                    }
                }

                if (group.ContainsKey(currentExtremum))
                {
                    continue;
                }

                group.Add(currentExtremum, currentGroup);

            }

            sharedExtrema = group.Where(p => p.Value.Count() >= count_threshold * XICs.Count())
                                 .Select(p=>p.Key).ToList();

        }

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

        public Dictionary<double, Tuple<double, double>> buildIndexedPeaks()
        {
            var indexedPeaks = new Dictionary<double, Tuple<double, double>>();

            foreach (var extremPoint in sharedExtrema.Where(p=>p.Type == ExtremumType.Maximum))
            {
                int index_Peak = sharedExtrema.IndexOf(extremPoint);
                int index_left = sharedExtrema
                    .Where(p => p.Type == ExtremumType.Minimum && p.RetentionTime < extremPoint.RetentionTime)
                    .Select(p => sharedExtrema.IndexOf(p))
                    .LastOrDefault();

                int index_right = sharedExtrema
                    .Where(p => p.Type == ExtremumType.Minimum && p.RetentionTime > extremPoint.RetentionTime)
                    .Select(p => sharedExtrema.IndexOf(p))
                    .FirstOrDefault();

                if(index_right == 0) index_right = sharedExtrema.Count() - 1;

                indexedPeaks.Add(extremPoint.RetentionTime, new Tuple<double, double> (sharedExtrema[index_left].RetentionTime, sharedExtrema[index_right].RetentionTime));
            }
            return indexedPeaks;
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
