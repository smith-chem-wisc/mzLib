using IO.MzML;
using NUnit.Framework;
using Proteomics;
using Chemistry;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace Test
{
    [TestFixture]
    public sealed class MyTests
    {
        #region Public Methods

        [Test]
        public void ModificationCollectionTest()
        {
            OldSchoolModification mod1 = new OldSchoolModification(10, "mass 10 modification");
            OldSchoolModification mod2 = new OldSchoolModification(100, "mass 100 modification");
            OldSchoolModification mod3 = new OldSchoolModification(1000, "mass 1000 modification");
            ModificationCollection a = new ModificationCollection(mod1, mod2, mod3, mod1);
            ModificationCollection b = new ModificationCollection(mod1, mod3, mod1, mod2);
            Assert.IsTrue(a.Equals(b));
            ModificationCollection c = new ModificationCollection(mod1);
            Assert.IsFalse(c.Equals(b));
        }

        [Test]
        public void MzmlReadTest()
        {
            //constructs three ms1 scans and three ms2 scans
            //if ms2 scan precusor reference is null, assumes most recent ms1 scan is correct reference
            // MzML\Mzml.cs @ lines 459-47

            MsDataScan[] scans = new MsDataScan[6];
            MsDataScan[] scans1 = new MsDataScan[4];

            double[] intensities0 = new double[] { 1 };
            double[] mz0 = new double[] { 50 };
            MzSpectrum massSpec0 = new MzSpectrum(mz0, intensities0, false);
            scans[0] = new MsDataScan(massSpec0, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec0.SumOfAllY, null, null, "1");
            scans1[0] = scans[0];

            double[] intensities1 = new double[] { 1 };
            double[] mz1 = new double[] { 50 };
            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            scans[1] = new MsDataScan(massSpec0, 2, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1");

            double[] intensities2 = new double[] { 1 };
            double[] mz2 = new double[] { 50 };
            MzSpectrum massSpec2 = new MzSpectrum(mz2, intensities2, false);
            scans[2] = new MsDataScan(massSpec2, 3, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec2.SumOfAllY, null, null, "1");

            //ms2
            double[] intensities3 = new double[] { 1 };
            double[] mz3 = new double[] { 30 };
            MzSpectrum massSpec3 = new MzSpectrum(mz3, intensities3, false);
            scans[3] = new MsDataScan(massSpec3, 4, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec3.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 3, null);
            scans1[1] = new MsDataScan(massSpec3, 2, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec3.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            //ms2
            double[] intensities4 = new double[] { 1 };
            double[] mz4 = new double[] { 30 };
            MzSpectrum massSpec4 = new MzSpectrum(mz4, intensities4, false);
            scans[4] = new MsDataScan(massSpec4, 5, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec4.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 3, null);
            scans1[2] = new MsDataScan(massSpec4, 3, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec4.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            //ms2
            double[] intensities5 = new double[] { 1 };
            double[] mz5 = new double[] { 30 };
            MzSpectrum massSpec5 = new MzSpectrum(mz5, intensities5, false);
            scans[5] = new MsDataScan(massSpec5, 6, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec5.SumOfAllY, null, null, "4", 50, null, null, 50, 1, DissociationType.CID, null, null);
            scans1[3] = new MsDataScan(massSpec5, 4, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec5.SumOfAllY, null, null, "4", 50, null, null, 50, 1, DissociationType.CID, null, null);

            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"), false);
            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"));

            FakeMsDataFile f1 = new FakeMsDataFile(scans1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f1, Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"), false);
            Mzml ok1 = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"));

            Assert.AreEqual(3, ok.GetAllScansList().ElementAt(5).OneBasedPrecursorScanNumber);
            Assert.AreEqual(1, ok1.GetAllScansList().ElementAt(3).OneBasedPrecursorScanNumber);
        }

        #endregion Public Methods
    }
}