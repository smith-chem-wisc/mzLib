using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using MzLibUtil;
using IO.Thermo;
using IO.MzML;
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
        public static void MergeBoxCarScans(MsDataFile file, SetOfBoxcarRanges[] boxcarRanges, string finalFilePath)
        {
            //SetOfBoxcarRanges[] boxcarRanges = FindBoxcars(file);
            List<SetOfScans> scans = SeparateScans(file);
            List<MsDataScan> mergedScans = MergeScans(scans, boxcarRanges);
            WriteMzmlFile(mergedScans, file, finalFilePath);
        }

        /// <summary>
        /// I adapted this method from Rob's Boxcar work
        /// loads a data file that has the extension .RAW as a ThermoDataFile
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns> ThermoDataFile
        public ThermoDataFile LoadThermoDataFile(string filepath)
        {
            string theExtension = Path.GetExtension(filepath).ToUpper();

            ThermoDataFile thermoDataFile = null;

            if (theExtension.Equals(".RAW"))
            {
                thermoDataFile = ThermoStaticData.LoadAllStaticData(filepath, null);
            }
            // can't load an mzml as a thermodatafile
            //else if (theExtension.Equals(".MZML"))
            //{
            //    thermoDataFile = Mzml.LoadAllStaticData(filepath);
            //}
            else
            {
                Console.WriteLine("Unrecognized format");
            }
            return thermoDataFile;
        }

        /// <summary>
        /// This method is from Rob's boxcar work
        /// loads a data file that has the extension .RAW or .MZML as an MsDataFile
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns> MsDataFile
        public MsDataFile LoadMsDataFile(string filepath)
        {
            string theExtension = Path.GetExtension(filepath).ToUpper();

            MsDataFile msDataFile = null;

            if (theExtension.Equals(".RAW"))
            {
                msDataFile = ThermoStaticData.LoadAllStaticData(filepath, null);
            }
            else if (theExtension.Equals(".MZML"))
            {
                msDataFile = Mzml.LoadAllStaticData(filepath);
            }
            else
            {
                Console.WriteLine("Unrecognized format");
            }
            return msDataFile;
        }


        // Finding boxcar ranges:

        /// <summary>
        /// Extracts boxcar m/z range information from the metadata (ThermoGlobalParams.InstrumentMethods),
        /// This information is stored in BoxcarRange objects.
        /// The BoxcarRange objects are stored in a SetOfBoxcarRanges, and each SetOfBoxcarRanges object has a number of BoxcarRange objects equal to the number of boxes specified in the scan setup.
        /// Depending on the number of boxcar scans (msx) taken, the SetOfBoxcarRanges[] may be a different length.
        /// (Currently this method is set up to return a SetOfBoxcarRanges[] of length 3, because there were 3 distinct boxcar scans).
        /// Removes overlap between boxcars to prepare for removal of edge effects when matched with the scan data.
        /// </summary>
        /// <param name="file"></param>
        /// <returns></returns> SetOfBoxcarRanges[] where each element in the array is a SetOfBoxcarRanges object.
        public static SetOfBoxcarRanges[] FindBoxcars(ThermoDataFile file)
        {
            string description = file.ThermoGlobalParams.InstrumentMethods[0];

            string se2 = "ScanEvent id=\"2\"";
            string se3 = "ScanEvent id=\"3\"";
            string se4 = "ScanEvent id=\"4\"";

            int index2 = description.IndexOf(se2);
            int index3 = description.IndexOf(se3);
            int index4 = description.IndexOf(se4);

            List<string> events = new List<string>();
            events.Add(description.Substring(index2, index3 - index2));
            events.Add(description.Substring(index3, index4 - index3));
            events.Add(description.Substring(index4));

            SetOfBoxcarRanges[] boxcars = new SetOfBoxcarRanges[3];

            for (int i = 0; i < boxcars.Count(); i++)
            {
                var e = events[i];
                SetOfBoxcarRanges set = new SetOfBoxcarRanges();
                string[] split = e.Split('[');
                string s = split[1].Split(']')[0];
                string[] points = s.Split('(');
                int len = points.Length;

                for (int j = 1; j < len; j++)
                {
                    string point = points[j];
                    string nums = point.Split(')')[0];
                    nums = nums.Trim();
                    string[] splitnums = nums.Split(',');
                    double startpoint = Convert.ToDouble(splitnums[0]);
                    double endpoint = Convert.ToDouble(splitnums[1]);
                    BoxcarRange bc = new BoxcarRange(startpoint, endpoint);
                    set.AddBoxcarRange(bc);
                }
                boxcars[i] = set;
            }
            SetOfBoxcarRanges[] finalBoxcars = RemoveOverlap(boxcars);
            return finalBoxcars;
        }

        /// <summary>
        /// Helper Method for FindBoxcars
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
                        boxcars[x].Insert(rangeA, y);
                        if (x < (numBoxcarScans - 1) && (y <= boxcars[x + 1].Count())) // if you aren't on the last column, insert rangeB into the next column
                        {
                            boxcars[x + 1].Insert(rangeB, y);
                        }
                        else if (x == (numBoxcarScans - 1) && (y != boxcars[0].Count())) // if you're on the last column, insert rangeB into the first column in the next row
                        {
                            boxcars[0].Insert(rangeB, y + 1);
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
        /// Helper method for RemoveEdgeEffects
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


        // Merging the scans and writing an mzml file with the merged scans
          
        /// <summary>
        /// Merges a List of SetOfScans objects into a List of MsDataScan objects by combining the scans taken with the boxcar method in each SetOfScans object into one MsDataScan.
        /// The new MsDataScans result from a process that matches their data with boxcar scan ranges and removes edge effects from the scans.
        /// </summary>
        /// <param name="toMatch"></param>
        /// <param name="boxcarTypes"></param>
        /// <returns></returns> MsDataFile containing the merged scans
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

        /// <summary>
        /// Writes an mzml file to the given filepath.
        /// Metadata for the file comes from the originalFile.
        /// </summary>
        /// <param name="mergedScans"></param>
        /// <param name="originalFile"></param>
        public static void WriteMzmlFile(List<MsDataScan> mergedScans, MsDataFile originalFile, string filepath)
        {
            MsDataScan[] combined = new MsDataScan[mergedScans.Count];
            //MsDataScan[] ms1scans = new MsDataScan[mergedScans.Count];

            if (mergedScans.Count == 0)
            {
                Console.WriteLine("Error! You have no merged scans.");
                return;
            }

            // put the merged boxcar scans into this array:
            for (int i = 0; i < mergedScans.Count; i++)
            {
                MsDataScan boxcarScan = mergedScans.ElementAt(i); ; 
                //MsDataScan ms1Scan = set.Ms1scans[0];
                
                combined[i] = boxcarScan;
                //ms1scans[i] = ms1Scan;
            }

            SourceFile sourceFile = originalFile.SourceFile;
            MsDataFile msDataFile = new MsDataFile(combined, sourceFile);
            //MsDataFile ms1DataFile = new MsDataFile(ms1scans, sourceFile);

            // Debugging: (change directory before using)
            //Console.WriteLine("scan5: ");
            //var scan5 = combined[5];
            //double[] xarray = scan5.MassSpectrum.XArray;
            //double[] yarray = scan5.MassSpectrum.YArray;
            //for (int i = 0; i < xarray.Length; i++)
            //{
            //    Console.WriteLine("point: " + xarray[i] + ", " + yarray[i]);
            //}

            //using (StreamWriter writer = new StreamWriter("C:\\Nicole\\test5.txt"))
            //{
            //    writer.WriteLine("m/z value " + "\t" + "intensity");
            //    int len = combined[5].MassSpectrum.XArray.Count();
            //    //Console.WriteLine(len);
            //    for (int i = 0; i < len; i++)
            //    {
            //        writer.WriteLine(combined[5].MassSpectrum.XArray[i] + "\t" + combined[5].MassSpectrum.YArray[i]);
            //    }
            //}
            
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, filepath + originalFile.SourceFile.FileName + "_output3_combined_boxcars.mzML", false);
            //MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(ms1DataFile, filepath + originalFile.SourceFile.FileName + "_output2_ms1s.mzML", false);
        }


        // NOT CURRENTLY USED:

        /// <summary>
        /// NOT CURRENTLY USED
        /// Returns a List<SetOfScans> where all SetOfScans objects contain the same number of MsDataScans that were taken with the boxcar method.
        /// Given a List<SetOfScans>, finds the average number of boxcar method scans for a set in that list, then
        /// rounds that average and discards any scans that do not have that number of boxcar method MsDataScans.
        /// </summary>
        /// <param name="list"></param>
        /// <returns></returns>
        //private List<SetOfScans> RemoveIncompleteSets(List<SetOfScans> list)
        //{
        //    int totalCount = 0;
        //    int numSets = 0;
        //    foreach (var setOfScans in list)
        //    {
        //        totalCount += setOfScans.Boxcars.Count;
        //        numSets++;
        //    }
        //    double averageCount = totalCount / numSets;
        //    double idealCount = Math.Round(averageCount);

        //    List<SetOfScans> cleanedSets = new List<SetOfScans>();
        //    foreach (var setOfScans in list)
        //    {
        //        if (setOfScans.Boxcars.Count.Equals((int)idealCount))
        //        {
        //            cleanedSets.Add(setOfScans);
        //        }
        //    }
        //    return cleanedSets;
        //}

        /// <summary>
        /// Helper Method for FindBoxcars
        /// Removes the edge effects from the boxcar ranges
        /// so that each boxcar shares a startpoint with the endpoint of the previous car
        /// and an endpoint with the startpoint of the next car.
        /// </summary>
        /// <param name="boxcars"></param>
        /// <returns></returns>
        //private SetOfBoxcarRanges[] RemoveEdgeEffects(SetOfBoxcarRanges[] boxcars)
        //{
        //    // order the boxcars by what the startpoint is
        //    //List<BoxcarRange> ordered = OrderBoxcars(boxcars);

        //    List<BoxcarRange> newBoxcars = new List<BoxcarRange>();

        //    if (ordered.Count > 1)
        //    {
        //        BoxcarRange previous = ordered[0];
        //        BoxcarRange current = ordered[0];
        //        BoxcarRange next = ordered[1];
        //        double prevmean = CalculateMean(current.GetEnd(), next.GetStart());

        //        newBoxcars.Add(new BoxcarRange(current.GetStart(), prevmean));

        //        // advance the current boxcar
        //        current = ordered[1];

        //        for (int i = 1; i < (ordered.Count - 1); i++)
        //        {
        //            // find the next boxcar
        //            next = ordered[i + 1];

        //            // calculate the mean of the endpoint of the current boxcar and the startpoint of the next boxcar
        //            double num1 = current.GetEnd();
        //            double num2 = next.GetStart();
        //            double mean = CalculateMean(num1, num2);

        //            // adds new boxcar objects to the list
        //            // where the startpoint of the new boxcar is the endpoint of the previous boxcar,
        //            // and the endpoint of the new boxcar is the mean of the endpoint of the current boxcar and the startpoint of the next boxcar.
        //            BoxcarRange bc = new BoxcarRange(prevmean, mean);
        //            //Console.WriteLine("start: " + bc.GetStart() + "end: " + bc.GetEnd());
        //            newBoxcars.Add(bc);
        //            previous = current;
        //            current = next;
        //            prevmean = mean;
        //        }
        //        return newBoxcars;
        //    }
        //    else
        //    {
        //        return null;
        //    }
        //}

        /// <summary>
        /// Helper method for RemoveEdgeEffects
        /// sorts the boxcars by startpoint in ascending order
        /// </summary>
        /// <param name="boxcars"></param>
        /// <returns></returns>
        //private List<BoxcarRange> OrderBoxcars(List<BoxcarRange> boxcars)
        //{
        //    List<BoxcarRange> newBoxcars = new List<BoxcarRange>();

        //    int count = boxcars.Count;
        //    BoxcarRange[] bcarray = new BoxcarRange[count];
        //    for (int i = 0; i < count; i++)
        //    {
        //        bcarray[i] = boxcars.ElementAt(i);
        //    }

        //    Array.Sort(bcarray, BoxcarRange.SortAscending());

        //    for (int i = 0; i < count; i++)
        //    {
        //        newBoxcars.Add(bcarray[i]);
        //    }

        //    //Console.WriteLine("post-ordering: ");
        //    //foreach (var b in newBoxcars)
        //    //{
        //    //    Console.WriteLine(b.GetStart() + "\t" + b.GetEnd());
        //    //}
        //    return newBoxcars;
        //}

        /// <summary>
        /// NOT CURRENTLY USED
        /// separates scans in a ThermoDataFile file into 3 categories (all boxcar scans), discards other scans
        /// categories: Full msx ms (a), Full msx ms (b), Full msx ms (c)
        /// For each group of 3 scans, creates a SetOfScans object and adds to an arraylist
        /// The method discards boxcar scans if they do not have a corresponding full ms1 scan.
        /// </summary>
        /// <param name="file"></param>
        /// <returns></returns>
        //private List<SetOfScans> GetBoxcarScans(ThermoDataFile file)
        //{
        //    // hold the SetOfScans objects
        //    List<SetOfScans> sorted = new List<SetOfScans>();

        //    // discard all ms2 scans
        //    var ms1scans = file.GetMS1Scans().ToList();

        //    int len = ms1scans.Count;

        //    // initialize variables
        //    scan0 = ms1scans.ElementAt(0);
        //    MsDataScan scan1 = scan0;
        //    MsDataScan scan2 = scan0;
        //    MsDataScan scan3 = scan0;

        //    scanFilter0 = "";
        //    string scanFilter1 = "";
        //    string scanFilter2 = "";
        //    string scanFilter3 = "";

        //    // loop through all scans, adding sets of 3 boxcar scans (msx, msx, msx) to the sorted arraylist
        //    for (int i = 0; i <= len; i++)
        //    {
        //        if (scanCounter < 4) // if you haven't just added a SetOfScans object to the sorted arraylist.
        //        {
        //            if (scanCounter == 0)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan0 = CheckMS(scan);

        //            }
        //            else if (scanCounter == 1)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan1 = CheckMSX(scan);
        //                scanFilter1 = scan1.ScanFilter;
        //            }
        //            else if (scanCounter == 2)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan2 = CheckMSX(scan);
        //                scanFilter2 = scan2.ScanFilter;
        //            }
        //            else if (scanCounter == 3)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan3 = CheckMSX(scan);
        //                scanFilter3 = scan3.ScanFilter;
        //            }
        //        }
        //        // if you have 3 boxcar "Full msx ms" scans, add them to a new setofscans object, add that object to the sorted list.
        //        else if (scanCounter >= 4)
        //        {
        //            if (scanFilter0.Contains("Full ms") && scanFilter1.Contains("Full msx ms") && scanFilter2.Contains("Full msx ms") && scanFilter3.Contains("Full msx ms"))
        //            {
        //                SetOfScans s = new SetOfScans(scan1, scan2, scan3);
        //                sorted.Add(s);
        //                //Console.WriteLine("Added new SetOfScans object: (1) " + scanFilter1 + " with id: " + scan1.NativeId + ", (2) " + scanFilter2 + " with id: " + scan2.NativeId + ", (3)" + scanFilter3 + " with id: " + scan3.NativeId + " to arraylist sorted");
        //            }
        //            if (i != len) // if there are scans left
        //            {
        //                scanCounter = 0;
        //                scan0 = CheckMS(ms1scans.ElementAt(i));
        //            }
        //        }
        //    }
        //    //Console.WriteLine("Length of sorted: " + sorted.Count);

        //    return sorted;
        //}

        /// <summary>
        /// NOT CURRENTLY USED
        /// Helper method for GetBoxcarScans
        /// Checks an MsDataScan to see if it is an msx (boxcar)
        /// If it is not a boxcar, checks to see if it is a regular ms1 scan and sets the scanCounter to 1, or to 0 if it is neither a boxcar or an ms1.
        /// If it is a boxcar, increments the scan count
        /// 
        /// </summary>
        /// <param name="scan"></param>
        /// <returns></returns>
        //private MsDataScan CheckMSX(MsDataScan scan)
        //{
        //    String scanFilter = scan.ScanFilter;
        //    scanCounter++;
        //    //Console.WriteLine("Scan Filter: " + scanFilter);

        //    // make sure the scan is a full msx ms (no scans are missing from the data)
        //    if (!scanFilter.Contains("Full msx ms"))
        //    {
        //        // check to see if the scan is a full ms, if so, record it and restart count from 1
        //        if (scanFilter.Contains("Full ms "))
        //        {
        //            scan0 = scan;
        //            scanFilter0 = scanFilter;
        //            scanCounter = 1;
        //            //Console.WriteLine("full ms, restarting count from 1");
        //        }
        //        // if the scan is neither a full ms or a full msx ms, do nothing and move on to the next scan
        //        else
        //        {
        //            scanCounter = 0;
        //            //Console.WriteLine("not a full msx, restarting count from 0");
        //        }
        //    }
        //    else
        //    {
        //        //Console.WriteLine("Success. Count incremented.");
        //    }
        //    return scan;
        //}

        /// <summary>
        /// NOT CURRENTLY USED
        /// Helper method for GetBoxcarScans
        /// Checks to see if the MsDataScan is a full ms1 scan. 
        /// If not, sets the scanCounter to 0.
        /// If it is, increments the count.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="ms1scans"></param>
        /// <returns></returns>
        //private MsDataScan CheckMS(MsDataScan scan)
        //{
        //    scan0 = scan;
        //    scanFilter0 = scan0.ScanFilter;
        //    scanCounter++;
        //    //Console.WriteLine("Scan Filter: " + scanFilter0);
        //    if (!scanFilter0.Contains("Full ms "))
        //    {
        //        scanCounter = 0;
        //        //Console.WriteLine("not a full ms, restarting count from 0");
        //    }
        //    else
        //    {
        //        //Console.WriteLine("Success. Count incremented.");
        //    }
        //    return scan0;
        //}

        /// <summary>
        /// NOT CURRENTLY USED
        /// separates scans in a ThermoDataFile file into 4 categories, discards other scans
        /// categories: Full ms, Full msx ms (a), Full msx ms (b), Full msx ms (c)
        /// For each group of 4 scans, creates a SetOfScans object and adds to an arraylist
        /// </summary>
        /// <param name="file"></param>
        //private List<SetOfScans> SeparateScans4(ThermoDataFile file)
        //{
        //    // hold the SetOfScans objects
        //    List<SetOfScans> sorted = new List<SetOfScans>();

        //    // discard all ms2 scans
        //    var ms1scans = file.GetMS1Scans().ToList();

        //    int len = ms1scans.Count;

        //    // initialization of variables
        //    scanCounter = 0;

        //    scan0 = ms1scans.ElementAt(0);
        //    MsDataScan scan1 = scan0;
        //    MsDataScan scan2 = scan0;
        //    MsDataScan scan3 = scan0;

        //    scanFilter0 = "";
        //    string scanFilter1 = "";
        //    string scanFilter2 = "";
        //    string scanFilter3 = "";

        //    // loop through all scans, adding full sets of 4 scans (ms, msx, msx, msx) to the sorted arraylist
        //    for (int i = 0; i < len; i++)
        //    {
        //        if (scanCounter < 4) // if you haven't just added a SetOfScans object to the sorted arraylist.
        //        {
        //            if (scanCounter == 0)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan0 = CheckMS(scan);

        //            }
        //            else if (scanCounter == 1)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan1 = CheckMSX(scan);
        //                scanFilter1 = scan1.ScanFilter;
        //            }
        //            else if (scanCounter == 2)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan2 = CheckMSX(scan);
        //                scanFilter2 = scan2.ScanFilter;
        //            }
        //            else if (scanCounter == 3)
        //            {
        //                MsDataScan scan = ms1scans.ElementAt(i);
        //                scan3 = CheckMSX(scan);
        //                scanFilter3 = scan3.ScanFilter;
        //            }
        //        }
        //        // if you have 4 scans, 1 "Full ms " and 3 "Full msx ms", add them to a new setofscans object, add that object to the sorted list.
        //        else if (scanCounter >= 4)
        //        {
        //            if (scanFilter0.Contains("Full ms") && scanFilter1.Contains("Full msx ms") && scanFilter2.Contains("Full msx ms") && scanFilter3.Contains("Full msx ms"))
        //            {
        //                SetOfScans s = new SetOfScans(scan0, scan1, scan2, scan3);
        //                sorted.Add(s);
        //                Console.WriteLine("Added new SetOfScans object: (0) " + scanFilter0 + " with id: " + scan0.NativeId + ", (1) " + scanFilter1 + " with id: " + scan1.NativeId + ", (2) " + scanFilter2 + " with id: " + scan2.NativeId + ", (3)" + scanFilter3 + " with id: " + scan3.NativeId + " to arraylist sorted");
        //            }
        //            scanCounter = 0;
        //            scan0 = CheckMS(ms1scans.ElementAt(i));
        //        }
        //    }
        //    Console.WriteLine("Length of sorted: " + sorted.Count);
        //    return sorted;
        //}      

        /// <summary>
        /// NOT CURRENTLY USED.
        /// MOSTLY FOR DEBUGGINGS
        /// Writes 4 files to the file structure, containing the x and y arrays for all the scans in the input set.
        /// probably want to put these in a folder right away.
        /// </summary>
        /// <param name="set"></param>
        //private void WriteFiles(SetOfScans set, int i)
        //{
        //    foreach (var scan in set)
        //    {
        //        if (scan.ScanFilter.Contains("Full ms "))
        //        {
        //            using (StreamWriter write = new StreamWriter("C:\\Nicole\\fullms.txt"))
        //            {
        //                int arrayLength = scan.MassSpectrum.XArray.Length;
        //                write.WriteLine(scan.ScanFilter);
        //                write.WriteLine("numRows = " + arrayLength);
        //                for (int j = 0; j < arrayLength; j++)
        //                {
        //                    write.WriteLine(scan.MassSpectrum.XArray[j] + "\t" + scan.MassSpectrum.YArray[j]);
        //                }
        //            }
        //        }
        //        else
        //        {
        //            using (StreamWriter write = new StreamWriter("C:\\Nicole\\boxcar" + i + ".txt"))
        //            {
        //                int arrayLength = scan.MassSpectrum.XArray.Length;
        //                write.WriteLine(scan.ScanFilter);
        //                write.WriteLine("numRows = " + arrayLength);
        //                for (int j = 0; j < arrayLength; j++)
        //                {
        //                    write.WriteLine(scan.MassSpectrum.XArray[j] + "\t" + scan.MassSpectrum.YArray[j]);
        //                }
        //            }
        //        }
        //    }
        //}
    }
}

