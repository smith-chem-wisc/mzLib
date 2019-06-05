using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using MzLibUtil;
using MassSpectrometry;

namespace BoxCar
{
    /// <summary>
    /// Merges the scans in a data file that were taken using the boxcar method.
    /// Writes an mzml file containing the merged scans.
    /// 
    /// Written by Nicole Frey, May-June 2019 for the Smith Group in the UW Madison chemistry department, with direction from Leah Schaffer.
    /// </summary>
    public class Program
    {
        /// <summary>
        /// Uses the boxcar algorithm to merge the boxcar scans in a file, writes a new mzml file to the directory.
        /// </summary>
        /// <param name="file"></param>
        /// <param name="boxcarRanges"></param>
        /// <param name="finalFilePath"></param>
        /// <returns></returns> mergedScans, a list of final merged boxcar scans.
        public static List<MsDataScan> MergeBoxCarScans(MsDataFile file, SetOfBoxcarRanges[] boxcarRanges, string finalFilePath)
        {
            //SetOfBoxcarRanges[] boxcarRanges = FindBoxcars(file);
            SetOfBoxcarRanges[] bcRanges = RemoveOverlap(boxcarRanges);
            List<SetOfScans> scans = SeparateScans(file);
            List<MsDataScan> mergedScans = MergeScans(scans, bcRanges);
            // WriteMzmlFile(mergedScans, file, finalFilePath);
            return mergedScans;
        }


        /// <summary>
        /// Removes the overlap from the boxcar ranges
        /// so that each boxcar shares a startpoint with the endpoint of the previous car
        /// and an endpoint with the startpoint of the next car.
        /// 
        /// Does this process by taking the mean between the boxcar ranges and setting the cars' start and endpoints to the means of their neighbors.
        /// Overlap removal is important because when these boxcar ranges are matched with the scan data, they will remove the edge effects from the scan.
        /// Edge effects: sometimes the scan picks up on peaks that it shouldn't near the beginning and end of that boxcar.
        /// 
        /// EXAMPLE:
        /// 
        /// boxcars[0]      boxcars[1]      boxcars[2] ... (if you took more kinds of boxcar scans they would be additional columns)
        /// (400, 426)      (424, 451)      (449, 476)
        /// (474, 501)      (499, 526)      (524, 551)      (these rows are the BoxcarRanges)
        ///     .               .               .
        ///     .               .               .
        ///     .               .               .
        /// (1124, 1151)    (1149, 1176)    (1174,1201) 
        /// 
        /// How the loop works:
        /// there are x columns and y rows.
        /// Each column is a SetOfBoxcarRanges (which is stored in an arraylist) 
        /// Each row contains some BoxcarRanges, and if you read across the rows or down the columns they increase.
        /// The nested for loops look at the rows first, and then the columns. 
        /// This way, each BoxcarRange that is looped through is greater than the last BoxcarRange and you can calculate the means.
        /// 
        /// So, once the function loops through, the entries are changed to:
        /// 
        /// boxcars[0]      boxcars[1]      boxcars[2] ... 
        /// (400, 425)      (425, 450)      (450, 475)
        /// (475, 500)      (500, 525)      (525, 550)      
        ///     .               .               .
        ///     .               .               .
        ///     .               .               .
        /// (1125, 1150)    (1150, 1175)    (1175,1201) - the last entry (and the first entry) doesn't change because there's no next entry to compute the mean
        /// 
        /// In reality, the intervals are not all the same, so it's slightly more complicated than the example
        /// 
        /// </summary>
        /// <param name="boxcars"></param>
        /// <returns></returns>
        public static SetOfBoxcarRanges[] RemoveOverlap(SetOfBoxcarRanges[] boxcars)
        {
            int numBoxcarScans = boxcars.Count();

            SetOfBoxcarRanges[] newBoxcars = new SetOfBoxcarRanges[numBoxcarScans];

            if (numBoxcarScans > 1)
            {
                BoxcarRange rangeA = boxcars[0].ElementAt(0);
                BoxcarRange rangeB = boxcars[1].ElementAt(0);
                // loop through all complete rows and columns of the "matrix":
                for (int y = 0; y < boxcars[0].Count(); y++) // loops through all rows including a possibly incomplete last row w/ nulls/empty cells in it
                {
                    for (int x = 0; x < numBoxcarScans; x++)
                    {
                        if (x < (numBoxcarScans - 1) && (y <= boxcars[x + 1].Count())) // if you aren't on the last column and there is an entry in the yth row of the next column, the rangeB is in the next column
                        {
                            rangeB = boxcars[x + 1].ElementAt(y);
                        }
                        else if (x == (numBoxcarScans - 1) && (y != (boxcars[0].Count() - 1))) // if you're on the last column and there is another row (even a partial row), rangeB is the first column in the next row
                        {
                            rangeB = boxcars[0].ElementAt(y + 1);
                        }
                        else // if you've reached the last entry
                        {
                            return boxcars;
                        }
                        // find the mean of the endpoint of rangeA and the startpoint of rangeB
                        double endA = rangeA.End;
                        double startB = rangeB.Start;
                        double mean = CalculateMean(endA, startB);
                        // change the endpoint of rangeA and the startpoint of rangeB to be that number
                        rangeA.End = mean;
                        rangeB.Start = mean;
                        // insert rangeA and rangeB
                        boxcars[x].ReplaceAtIndex(rangeA, y);
                        if (x < (numBoxcarScans - 1) && (y <= boxcars[x + 1].Count())) // if you aren't on the last column, insert rangeB into the next column
                        {
                            boxcars[x + 1].ReplaceAtIndex(rangeB, y);
                        }
                        else if (x == (numBoxcarScans - 1) && (y != boxcars[0].Count())) // if you're on the last column, insert rangeB into the first column in the next row
                        {
                            boxcars[0].ReplaceAtIndex(rangeB, y + 1);
                        }
                        rangeA = rangeB;
                    }
                }
                return boxcars;
            }
            else
            {
                return boxcars;
            }
        }

        /// <summary>
        /// Helper method for RemoveOverlap
        /// given 2 numbers, returns the mean
        /// </summary>
        /// <param name="num1"></param>
        /// <param name="num2"></param>
        /// <returns></returns> double mean
        public static double CalculateMean(double num1, double num2)
        {
            return ((num1 + num2) / 2);
        }


        // Grouping the scans in the file

        /// <summary>
        /// Given a file, reads through the scans and groups them into SetOfScans objects, each containing one or more full ms1 scans and the corresponding boxcar scans.
        /// Discards all ms2 scans (although the SetOfScans object also has a list to hold these, if you wanted to include them).
        /// This method includes many comments which can be uncommented for debugging to write information to the console.
        /// </summary>
        /// <param name="file"></param>
        /// <returns></returns> list of SetOfScans objects, each SetOfScans contains one full ms1 scan and its corresponding boxcar scans.
        public static List<SetOfScans> SeparateScans(MsDataFile file)
        {
            List<SetOfScans> sorted = new List<SetOfScans>();

            // Get only the ms1 and boxcar scans (discard ms2):
            var ms1scans = file.GetMS1Scans().ToList();
            //Console.WriteLine("Number of ms1scans (full ms and boxcar method): " + ms1scans.Count);

            // Initialize variables:

            MsDataScan current = ms1scans.ElementAt(0);
            string currentFilter = current.ScanFilter;
            //Console.WriteLine("first scan filter: " + currentFilter);

            MsDataScan next = ms1scans.ElementAt(1);
            string nextFilter = next.ScanFilter;
            //Console.WriteLine("next scan filter: " + nextFilter);

            SetOfScans set = new SetOfScans();

            // Loop through to group the scans. This should create as many SetOfScans objects as there are boxes (specified in the setup of the scan).
            for (int i = 0; i < ms1scans.Count; i++)
            {
                if (current.ScanFilter.Contains("Full ms ")) // the space is important! Otherwise this will be recorded as a boxcar scan instead of an ms1
                {
                    set.AddToMs1Scans(current);
                   // Console.WriteLine("Found full ms, added to set");
                }
                // if this boxcar scan has boxcar scans after it
                else if (current.ScanFilter.Contains("Full msx ms") && next.ScanFilter.Contains("Full msx ms"))
                {
                    set.AddToBoxcarScans(current);
                    //Console.WriteLine("Found msx ms, added to set");
                }
                // if this boxcar scan is the last boxcar scan in the set, start a new set
                else if (current.ScanFilter.Contains("Full msx ms") && !next.ScanFilter.Contains("Full msx ms"))
                {
                    set.AddToBoxcarScans(current);
                    //Console.WriteLine("Found last msx ms in set, added to set");
                    sorted.Add(set);
                    //Console.WriteLine("Added set to Sorted");
                    set = new SetOfScans();
                }
                current = next;
                currentFilter = current.ScanFilter;
                if ((i+1) < ms1scans.Count) // if there is a next MsDataScan in the file
                {
                    next = ms1scans.ElementAt(i + 1);
                    nextFilter = next.ScanFilter;
                }
            }
            return sorted;
        }


        // Merging the scans
          
        /// <summary>
        /// Merges a List of SetOfScans objects into a List of MsDataScan objects by combining the scans taken with the boxcar method in each SetOfScans object into one MsDataScan.
        /// The new MsDataScans result from a process that matches their data with boxcar scan ranges and removes edge effects from the scans.
        /// </summary>
        /// <param name="toMatch"></param>
        /// <param name="boxcarTypes"></param>
        /// <returns></returns> List of MsDatascans containing the merged scans
        public static List<MsDataScan> MergeScans(List<SetOfScans> toMatch, SetOfBoxcarRanges[] boxcarTypes)
        {
            List<MsDataScan> combinedScans = new List<MsDataScan>();

            foreach (var set in toMatch)
            {
                // combine all boxcar method scans into one scan
                MsDataScan combinedScan = CombineScans(set, boxcarTypes);

                // Code for returning a List<SetOfScans> instead of a List<MsDataScan>, so that the original ms1 (and ms2 if they are in the original toMatch List) stay with the combinedScan.
                // Not sure how well this works, may need to fix it/ test it to make sure.
                // since the new SetOfScans object we will put in combinedScans needs a List of boxcar-method scans, add newScan to a list (this will be a list of length one)
                //List<MsDataScan> newList = new List<MsDataScan>();
                //newList.Add(newScan);
                // create a new SetOfScans to hold the previous ms1 and ms2 information, as well as the combined boxcar scan.
                //SetOfScans newSet = new SetOfScans(set.Ms1scans, newList, set.Ms2scans);

                combinedScans.Add(combinedScan);
            }
            return combinedScans;
        }

        /// <summary>
        /// Helper method for MergeScans
        /// Combines the 3 scans in a set into one MsDataScan.
        /// </summary>
        /// <param name="set"></param>
        /// <param name="boxcarRanges"></param>
        /// <returns></returns> MsDataScan combined scans
        public static MsDataScan CombineScans(SetOfScans set, SetOfBoxcarRanges[] boxcarTypes)
        {
            List<double[]> bestPoints = new List<double[]>();
            var count = 1;

            foreach (var boxcarRanges in boxcarTypes)
            {
                foreach (var boxcarRange in boxcarRanges)
                {
                    BestScan bestScan = FindBoxcarScan(set, boxcarRange, count);

                    // find the m/z values and intensities in that scan's x and y arrays between the first and last index, add to the lists
                    for (int i = bestScan.FirstIndex; i <= bestScan.LastIndex; i++)
                    {
                        double x = bestScan.Scan.MassSpectrum.XArray[i];
                        double y = bestScan.Scan.MassSpectrum.YArray[i];
                        bestPoints.Add(new double[] { x, y });
                    }
                }
            }

            // recast bestPoints (List<Double[]>) into double[,]
            double[,] mzIntensities = new double[bestPoints[0].Length, bestPoints.Count];

            for (int i = 0; i < bestPoints[0].Length; i++)
            {
                for (int j = 0; j < bestPoints.Count; j++)
                {
                    mzIntensities[i, j] = bestPoints[j][i];
                }
            }

            // create and return the new MsDataScan
            MzSpectrum newMassSpectrum = new MzSpectrum(mzIntensities);
            MsDataScan newScan = new MsDataScan(newMassSpectrum, count, 1, true, set.BoxcarScans[0].Polarity, set.BoxcarScans[0].RetentionTime, set.BoxcarScans[0].ScanWindowRange, set.BoxcarScans[0].ScanFilter, set.BoxcarScans[0].MzAnalyzer, set.BoxcarScans[0].TotalIonCurrent, set.BoxcarScans[0].InjectionTime, set.BoxcarScans[0].NoiseData, set.BoxcarScans[0].NativeId);
            count++;
            return newScan;
        }

        /// <summary>
        /// Helper method for CombineScans
        /// Identifies which scan in the set has the most points in the boxcar range
        /// </summary>
        /// <param name="set"></param>
        /// <param name="boxcar"></param> BestScan object with first index, last index, total intensity, and scan
        public static BestScan FindBoxcarScan(SetOfScans set, BoxcarRange range, int count)
        {
            double bestIntensity = 0;
            MsDataScan bestScan = set.BoxcarScans[0];
            int bestFirstIndex = 0;
            int bestLastIndex = 0;
            int bestLen = 0;

            List<MsDataScan> boxcarScans = set.BoxcarScans;
            // Check each scan to find the one that corresponds best to the boxcar scan range
            foreach (var scan in boxcarScans)
            {
                BestScan bestrange = FindScanRange(scan, range);

                // the scan is better than the previous scan if:
                // it has more m/z values in the range, and if
                // the total intensity is higher

                int firstIndex = bestrange.FirstIndex;
                int lastIndex = bestrange.LastIndex;
                int currentLen = lastIndex - firstIndex;

                if (currentLen > bestLen)
                {
                    if (bestrange.TotalIntensity < bestIntensity)
                    {
                        Console.WriteLine("ACK! scan#: " + count + " has a lot of very low intensity scans and may have been incorrectly counted in a boxcar range!");
                    }

                    // if the scan has more m/z points in the boxcar range than the previously best scan, update variables.
                    bestLen = currentLen;
                    bestLastIndex = lastIndex;
                    bestFirstIndex = firstIndex;
                    bestScan = scan;
                }
            }
            return new BestScan(bestFirstIndex, bestLastIndex, bestIntensity, bestScan);
        }

        /// <summary>
        /// Helper method for FindBoxcarScan
        /// Finds the first and last index in the given scan where the m/z value in the scan is between the startpoint and endpoint of the given boxcar.
        /// The m/z value at the first index will be > the startpoint of the boxcar.
        /// The m/z value at the last index will be >= the endpoint of the boxcar.
        /// </summary>
        /// <param name="scan"></param>
        /// <param name="boxcar"></param>
        /// <returns></returns> BestScan object with a first index, last index, and total intensity of the scan.
        public static BestScan FindScanRange(MsDataScan scan, BoxcarRange range)
        {
            double[] xArray = scan.MassSpectrum.XArray;
            double[] yArray = scan.MassSpectrum.YArray;

            int len = xArray.Count();

            // find the first index where the m/z values of the scan are in the range of the boxcar
            int firstIndex = 0;
            double firstValue = xArray[firstIndex];
            bool nextIndexExists = true;
            while ((firstValue < range.Start) && nextIndexExists)
            {
                firstIndex++;
                firstValue = xArray[firstIndex];
                if (!(firstIndex + 1 < len))
                {
                    nextIndexExists = false;
                }
            }

            // find the last index where the m/z values of the scan are in the range of the boxcar
            int lastIndex = firstIndex;
            double totalIntensity = yArray[lastIndex];

            if (nextIndexExists)
            {
                int nextIndex = lastIndex + 1;
                double nextValue = xArray[nextIndex];

                while (nextValue <= range.End && nextIndexExists)
                {
                    lastIndex++;
                    nextIndex++;
                    totalIntensity += yArray[lastIndex];
                    if (nextIndex < len)
                    {
                        nextValue = xArray[nextIndex];
                    }
                    else
                    {
                        nextIndexExists = false;
                    }
                }
            }
            return new BestScan(firstIndex, lastIndex, totalIntensity, scan);
        }
    }
}

