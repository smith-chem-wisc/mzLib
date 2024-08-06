using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics.AminoAcidPolymer;
using System.IO;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using TopDownProteomics;
using FlashLFQ.ResultsReading;

namespace TestFlashLFQ
{
    internal class TestIdentificationAdapter
    {
        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv")]
        public void TestAddProteinGroupInfoCorrect(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            var allResults = file.ToList();

            List<Identification> identifications = new List<Identification>();
            identifications = IdentificationAdapter.MakeIdentifications(file);

            // list should contain five elements
            Assert.That(identifications.Count, Is.EqualTo(5));
            // one protein associated with given results, list should only contain this one element 
            Assert.That(identifications[0].ProteinGroups.Count, Is.EqualTo(1));
            // two proteins associated with given results, list should contain two elements
            Assert.That(identifications[2].ProteinGroups.Count, Is.EqualTo(2));
            
            Identification identification1= identifications[0];
            Assert.That(identification1.BaseSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(identification1.ModifiedSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(identification1.Ms2RetentionTimeInMinutes, Is.EqualTo(1.9398));
            Assert.That(identification1.MonoisotopicMass, Is.EqualTo(669.4173));
            Assert.That(identification1.PrecursorChargeState, Is.EqualTo(2));

            HashSet<ProteinGroup> proteinGroups = identification1.ProteinGroups;
            ProteinGroup proteinGroup1 = proteinGroups.First();
            Assert.That(proteinGroup1.ProteinGroupName, Is.EqualTo("P16403"));
            Assert.That(proteinGroup1.GeneName, Is.EqualTo("H12"));
            Assert.That(proteinGroup1.Organism, Is.EqualTo("HUMAN"));

            Identification identification5 = identifications[4];
            Assert.That(identification5.BaseSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(identification5.ModifiedSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(identification5.Ms2RetentionTimeInMinutes, Is.EqualTo(19.114));
            Assert.That(identification5.MonoisotopicMass, Is.EqualTo(724.398));
            Assert.That(identification5.PrecursorChargeState, Is.EqualTo(2));
        }
    }
}
