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
    public class TestBoxCar // does this need to be a sealed class?
    {
        static string FilepathMZML;
        static MsDataFile file;
        static SetOfBoxcarRanges[] boxcarRanges;

        [SetUp]
        public static void SetUp()
        {
            FilepathMZML = Path.Combine(TestContext.CurrentContext.TestDirectory, "20170802_QEp1_FlMe_SA_BOX0_SILAC_BoxCar_SLICED.mzML");
            file = Mzml.LoadAllStaticData(FilepathMZML, null);

            // Create boxcars

            ArrayList toAddA = new ArrayList();
            toAddA.Add(new BoxcarRange(400, 416.3));
            toAddA.Add(new BoxcarRange(441.2, 454.2));
            toAddA.Add(new BoxcarRange(476.3, 488.8));
            toAddA.Add(new BoxcarRange(510.3, 523.3));
            toAddA.Add(new BoxcarRange(545, 557.8));
            toAddA.Add(new BoxcarRange(580.8, 594));
            toAddA.Add(new BoxcarRange(618.4, 633));
            toAddA.Add(new BoxcarRange(660.3, 676.4));
            toAddA.Add(new BoxcarRange(708.3, 726.3));
            toAddA.Add(new BoxcarRange(764.4, 788.4));
            toAddA.Add(new BoxcarRange(837.9, 868.8));
            toAddA.Add(new BoxcarRange(945, 999));

            ArrayList toAddB = new ArrayList();
            toAddB.Add(new BoxcarRange(415.3, 429.7));
            toAddB.Add(new BoxcarRange(453.2, 465.9));
            toAddB.Add(new BoxcarRange(487.8, 499.9));
            toAddB.Add(new BoxcarRange(522.3, 534.8));
            toAddB.Add(new BoxcarRange(556.8, 569.6));
            toAddB.Add(new BoxcarRange(593, 606.6));
            toAddB.Add(new BoxcarRange(632, 646.8));
            toAddB.Add(new BoxcarRange(675.4, 692.3));
            toAddB.Add(new BoxcarRange(725.3, 745));
            toAddB.Add(new BoxcarRange(787.4, 812.4));
            toAddB.Add(new BoxcarRange(867.8, 903.5));
            toAddB.Add(new BoxcarRange(998, 1071.1));

            ArrayList toAddC = new ArrayList();
            toAddC.Add(new BoxcarRange(428.7, 442.2));
            toAddC.Add(new BoxcarRange(464.9, 477.3));
            toAddC.Add(new BoxcarRange(498.9, 511.3));
            toAddC.Add(new BoxcarRange(533.8, 546));
            toAddC.Add(new BoxcarRange(568.6, 581.8));
            toAddC.Add(new BoxcarRange(605.6, 619.4));
            toAddC.Add(new BoxcarRange(645.8, 661.3));
            toAddC.Add(new BoxcarRange(691.3, 709.3));
            toAddC.Add(new BoxcarRange(744, 765.4));
            toAddC.Add(new BoxcarRange(811.4, 838.9));
            toAddC.Add(new BoxcarRange(902.5, 946));
            toAddC.Add(new BoxcarRange(1070.1, 1201));

            boxcarRanges = new SetOfBoxcarRanges[3] { new SetOfBoxcarRanges(toAddA), new SetOfBoxcarRanges(toAddB), new SetOfBoxcarRanges(toAddC) };
        }

        // Test methods:

        [Test]
        public static void TestCalculateMean()
        {
            Assert.AreEqual(511.5, Program.CalculateMean(509, 514));
        }

        [Test]
        public static void TestRemoveOverlap()
        {
            SetOfBoxcarRanges[] boxcarRangesOverlapRemoved = Program.RemoveOverlap(boxcarRanges);
            Assert.AreEqual(453.7, boxcarRangesOverlapRemoved[1].ElementAt(1).Start);
            Assert.IsTrue(453.700001 > boxcarRangesOverlapRemoved[1].ElementAt(1).Start && 453.699999 < boxcarRangesOverlapRemoved[1].ElementAt(1).Start);

            // make sure it works even if the SetOfBoxcarRanges[] is of length one.
            SetOfBoxcarRanges[] boxcarRangesMini = new SetOfBoxcarRanges[1];
            boxcarRangesMini[0] = boxcarRanges[0];
            boxcarRangesMini = Program.RemoveOverlap(boxcarRangesMini);
            Assert.AreEqual(boxcarRanges[0], boxcarRangesMini[0]);

            // make sure that the means are right: the new boxcarranges should have the same start and endpoints
            Assert.AreEqual(boxcarRangesOverlapRemoved[0].ElementAt(0).End, boxcarRangesOverlapRemoved[1].ElementAt(0).Start);
        }

        [Test]
        public static void TestMergeBoxcarScans()
        {
            // Sort the boxcar scans into categories

            boxcarRanges = Program.RemoveOverlap(boxcarRanges);

            List<SetOfScans> scans = Program.SeparateScans(file);
            Assert.AreNotEqual(null, scans);
            Assert.AreNotEqual(0, scans.Count);
            Assert.IsTrue((400.861480 > scans.ElementAt(11).BoxcarScans[0].MassSpectrum.XArray[8] - 0.00001) && (400.861480 < scans.ElementAt(11).BoxcarScans[0].MassSpectrum.XArray[8] + 0.00001));


            // Create mzml
            List<MsDataScan> mergedScans = Program.MergeScans(scans, boxcarRanges);
            Assert.AreNotEqual(0, mergedScans.Count);
            Assert.IsTrue((403.23065185 > mergedScans.ElementAt(5).MassSpectrum.XArray[22] - 0.00001) && (403.23065185 < mergedScans.ElementAt(5).MassSpectrum.XArray[22] + 0.00001));
            Assert.AreEqual(275908.875, Math.Round(mergedScans.ElementAt(5).MassSpectrum.YArray[37], 3));
            //WriteMzmlFile(mergedScans, file, FinalFilePath);

            int len = mergedScans.Count;
            for (int i = 0; i < len; i++)
            {
                Assert.AreEqual(mergedScans.ElementAt(i).MassSpectrum.XArray, Program.MergeBoxCarScans(file, boxcarRanges, null).ElementAt(i).MassSpectrum.XArray);
            }
        }

        // Tests for the SetOfScans class:

        [Test]
        public static void TestSetOfScansConstructor1()
        {
            List<MsDataScan> list = file.GetAllScansList();
            SetOfScans set = new SetOfScans(list, list, list);
            Assert.AreEqual(list, set.BoxcarScans);
            Assert.AreEqual(list, set.Ms1scans);
            Assert.AreEqual(list, set.Ms2scans);
        }

        [Test]
        public static void TestSetOfScansConstructor2()
        {
            List<MsDataScan> list = file.GetAllScansList();
            SetOfScans set = new SetOfScans(list, list);
            Assert.AreEqual(list, set.BoxcarScans);
            Assert.AreEqual(list, set.Ms1scans);
        }

        [Test]
        public static void TestAddToBoxcarScans()
        {
            MsDataScan boxcarscanExample = file.GetOneBasedScan(2);
            SetOfScans setOfScansExample = new SetOfScans();
            setOfScansExample.AddToBoxcarScans(boxcarscanExample);
            Assert.AreEqual(boxcarscanExample, setOfScansExample.BoxcarScans.ElementAt(0));
        }

        [Test]
        public void TestAddToMs1Scans()
        {
            MsDataScan ms1scanExample = file.GetOneBasedScan(1);
            SetOfScans setOfScansExample = new SetOfScans();
            setOfScansExample.AddToMs1Scans(ms1scanExample);
            Assert.AreEqual(ms1scanExample, setOfScansExample.Ms1scans.ElementAt(0));
        }

        [Test]
        public static void TestAddToMs2Scans()
        {
            MsDataScan ms2scanExample = file.GetOneBasedScan(5);
            SetOfScans setOfScansExample = new SetOfScans();
            setOfScansExample.AddToMs2Scans(ms2scanExample);
            Assert.AreEqual(ms2scanExample, setOfScansExample.Ms2scans.ElementAt(0));
        }


        // Tests for the SetOfBoxcarRanges Class:

        [Test]
        public static void TestSetOfBoxcarRangesConstructor()
        {
            ArrayList a = new ArrayList();
            a.Add(new BoxcarRange(0, 90));
            a.Add(new BoxcarRange(120, 333.3));

            SetOfBoxcarRanges set = new SetOfBoxcarRanges(a);
            Assert.AreEqual(a, set.Set);
        }

        [Test]
        public static void TestSetOfboxcarRangesGenericConstructor()
        {
            SetOfBoxcarRanges set = new SetOfBoxcarRanges();
            Assert.AreEqual(new ArrayList(), set.Set);
        }

        [Test]
        public static void TestBoxcarRangeMethods()
        {
            SetOfBoxcarRanges set = new SetOfBoxcarRanges();
            BoxcarRange br1 = new BoxcarRange(0, 90);
            BoxcarRange br2 = new BoxcarRange(89, 121);
            BoxcarRange br3 = new BoxcarRange(120, 333.3);
            set.AddBoxcarRange(br1);
            set.AddBoxcarRange(br3);
            Assert.AreEqual(br3, set.ElementAt(1));
            set.ReplaceAtIndex(br2, 1);
            Assert.AreEqual(br1, set.ElementAt(0));
            Assert.AreEqual(br2, set.ElementAt(1));
            Assert.AreEqual(2, set.Count());
        }

        // Tests for BoxcarRange class

        [Test]
        public static void TestBoxcarRangeConstruction()
        {
            BoxcarRange br = new BoxcarRange();
            BoxcarRange br1 = new BoxcarRange(401, 418.5);
            br.Start = 401;
            br.End = 418.5;
            Assert.AreEqual(br.Start, br1.Start);
            Assert.AreEqual(br.End, br1.End);

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
