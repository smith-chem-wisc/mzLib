using System.Collections;
using IO.Thermo;
using ManagedThermoHelperLayer;
using System.IO;
using IO.MzML;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using BoxCar;



/// <summary>
/// by Nicole Frey
/// </summary>
namespace Test
{
    [TestFixture]
    public class TestBoxCar // what does sealed mean??
    {
        string FilepathMZML = Path.Combine(TestContext.CurrentContext.TestDirectory, "20170802_QEp1_FlMe_SA_BOX0_SILAC_BoxCar_SLICED.mzML");

        [Test]
        public void TestCalculateMean()
        {
            Assert.AreEqual(511.5, Program.CalculateMean(509, 514));
            // maybe make sure this works with negative numbers OR make some sort of exception if negative numbers are inputted?
        }

        [Test]
        public void TestMergeBoxCarScans()
        {
            // Import data

            MsDataFile file = Mzml.LoadAllStaticData(FilepathMZML, null);

            ArrayList toAdd = new ArrayList();
            toAdd.Add(new BoxcarRange(400, 416.3));
            toAdd.Add(new BoxcarRange(441.2, 454.2));
            toAdd.Add(new BoxcarRange(476.3, 488.8));
            toAdd.Add(new BoxcarRange(510.3, 523.3));
            toAdd.Add(new BoxcarRange(545, 557.8));
            toAdd.Add(new BoxcarRange(580.8, 594));
            toAdd.Add(new BoxcarRange(618.4, 633));
            toAdd.Add(new BoxcarRange(660.3, 676.4));
            toAdd.Add(new BoxcarRange(708.3, 726.3));
            toAdd.Add(new BoxcarRange(764.4, 788.4));
            toAdd.Add(new BoxcarRange(837.9, 868.8));
            toAdd.Add(new BoxcarRange(945, 999));

            ArrayList toAddB = new ArrayList();
            toAdd.Add(new BoxcarRange(415.3, 429.7));
            toAdd.Add(new BoxcarRange(453.2, 465.9));
            toAdd.Add(new BoxcarRange(487.8, 499.9));
            toAdd.Add(new BoxcarRange(522.3, 534.8));
            toAdd.Add(new BoxcarRange(556.8, 569.6));
            toAdd.Add(new BoxcarRange(593, 606.6));
            toAdd.Add(new BoxcarRange(632, 646.8));
            toAdd.Add(new BoxcarRange(675.4, 692.3));
            toAdd.Add(new BoxcarRange(725.3, 745));
            toAdd.Add(new BoxcarRange(787.4, 812.4));
            toAdd.Add(new BoxcarRange(867.8, 903.5));
            toAdd.Add(new BoxcarRange(998, 1071.1));

            ArrayList toAddC = new ArrayList();
            toAdd.Add(new BoxcarRange(428.7, 442.2));
            toAdd.Add(new BoxcarRange(464.9, 477.3));
            toAdd.Add(new BoxcarRange(498.9, 511.3));
            toAdd.Add(new BoxcarRange(533.8, 546));
            toAdd.Add(new BoxcarRange(568.6, 581.8));
            toAdd.Add(new BoxcarRange(605.6, 619.4));
            toAdd.Add(new BoxcarRange(645.8, 661.3));
            toAdd.Add(new BoxcarRange(691.3, 709.3));
            toAdd.Add(new BoxcarRange(744, 765.4));
            toAdd.Add(new BoxcarRange(811.4, 838.9));
            toAdd.Add(new BoxcarRange(902.5, 946));
            toAdd.Add(new BoxcarRange(1070.1, 1201));

            SetOfBoxcarRanges[] boxcarRanges = new SetOfBoxcarRanges[3] { new SetOfBoxcarRanges(toAdd), new SetOfBoxcarRanges(toAddB), new SetOfBoxcarRanges(toAddC) } ;



            // Sort the boxcar scans into categories
            List<SetOfScans> scans = Program.SeparateScans(file);
            Assert.AreNotEqual(null, scans);
            Assert.AreNotEqual(0, scans.Count);
            Assert.IsTrue((400.861480 > scans.ElementAt(11).BoxcarScans[0].MassSpectrum.XArray[8] - 0.00001) && (400.861480 < scans.ElementAt(11).BoxcarScans[0].MassSpectrum.XArray[8] + 0.00001));


            // Create mzml
            List <MsDataScan> mergedScans = Program.MergeScans(scans, boxcarRanges);
            Assert.AreNotEqual(0, mergedScans.Count);
            Assert.IsTrue((403.23065185 > mergedScans.ElementAt(5).MassSpectrum.XArray[22] - 0.00001) && (403.23065185 < mergedScans.ElementAt(5).MassSpectrum.XArray[22] + 0.00001));
            //Assert.IsTrue(()) // 5 39 275908.875
            Assert.AreEqual(275908.875, Math.Round(mergedScans.ElementAt(5).MassSpectrum.YArray[37], 3));
            //WriteMzmlFile(mergedScans, file, FinalFilePath);

            //Assert.AreEqual(3, boxcarRanges.Count());
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

        //[Test]
        //public void TestFindBoxcars()
        //{
        //    // not sure if i need to test this one because i don't think it would be used like this in our data
        //}

        //[Test]
        //public void TestRemoveOverlap()
        //{
        //    // make sure that it returns the input SetOfBoxcarRanges[] if there is 0 or 1 in the array (numBoxcarScans <= 1)
        //    // make sure that the means are right:
        //    // the new boxcarranges should have the same start and endpoints
        //    // make sure it works even if the SetOfBoxcarRanges objects have different lengths (doubt this would be an issue, maybe test it anyway)
        //}

        //[Test]
        //public void TestSeparateScans()
        //{
        //    // make sure there are as many SetOfScans objects in the list as the scanheader # of boxes
        //}

        //[Test]
        //public void TestMergeScans()
        //{
        //    // returned list should have the same nubmer of elements as the input set of scans? maybe unless some got cut out. check that too
        //}

        //[Test]
        //public void TestCombineScans()
        //{
        //    // returned MsDataScan MsDataScan should have the same range as the boxcar scans, should be a combination
        //}

        //[Test]
        //public void TestFindBoxcarScan()
        //{
        //    // 
        //}

        //[Test]
        //public void TestFindScanRange()
        //{

        //}

        //[Test]
        //public void TestWriteMzmlFile()
        //{

        //}

        // Leah's data - manually finding boxcars because they aren't in the metadata:

        //ArrayList msxA = new ArrayList();
        //ArrayList msxB = new ArrayList();

        //msxA.Add(new BoxcarRange(500, 503.5));
        //msxA.Add(new BoxcarRange(505, 509));
        //msxA.Add(new BoxcarRange(510, 515));
        //msxA.Add(new BoxcarRange(518, 524));
        //msxA.Add(new BoxcarRange(528, 536));
        //msxA.Add(new BoxcarRange(544, 555));
        //msxA.Add(new BoxcarRange(567, 586));
        //msxA.Add(new BoxcarRange(607, 637));
        //msxA.Add(new BoxcarRange(673.4, 723));
        //msxA.Add(new BoxcarRange(782, 859));
        //msxA.Add(new BoxcarRange(952, 1071));
        //msxA.Add(new BoxcarRange(1214, 1391));

        //msxB.Add(new BoxcarRange(502.5, 506));
        //msxB.Add(new BoxcarRange(507.5, 511));
        //msxB.Add(new BoxcarRange(514, 519));
        //msxB.Add(new BoxcarRange(523, 529));
        //msxB.Add(new BoxcarRange(535, 545));
        //msxB.Add(new BoxcarRange(554, 568));
        //msxB.Add(new BoxcarRange(585, 608));
        //msxB.Add(new BoxcarRange(636, 674.5));
        //msxB.Add(new BoxcarRange(722, 783));
        //msxB.Add(new BoxcarRange(858, 954));
        //msxB.Add(new BoxcarRange(1070, 1215));
        //msxB.Add(new BoxcarRange(1390, 1600));

        //SetOfBoxcarRanges setA = new SetOfBoxcarRanges(msxA);
        //SetOfBoxcarRanges setB = new SetOfBoxcarRanges(msxB);

        //SetOfBoxcarRanges[] boxcarRanges = new SetOfBoxcarRanges[] { setA, setB };
        



  



   


    }
}
