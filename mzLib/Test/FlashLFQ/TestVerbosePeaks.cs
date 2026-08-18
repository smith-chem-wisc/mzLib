using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

namespace Test.FlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestVerbosePeaks
    {
        private static SpectraFileInfo SlicedMzml => new SpectraFileInfo(
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", "sliced-mzml.mzml"),
            "a", 0, 0, 0);

        private static FlashLfqResults Quantify(bool writeVerbosePeaks, out SpectraFileInfo mzml)
        {
            mzml = SlicedMzml;
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2,
                new List<ProteinGroup> { pg });

            var parameters = new FlashLfqParameters
            {
                MaxThreads = 1,
                Silent = true,
                WriteVerbosePeaks = writeVerbosePeaks
            };

            return new FlashLfqEngine(parameters, new List<Identification> { id }).Run();
        }

        /// <summary>
        /// Off by default, and off means the envelopes are the ordinary ones: nothing retains the isotope peaks, and
        /// the peaks file is written exactly as it was before the option existed.
        /// </summary>
        [Test]
        public static void TestVerbosePeaksAreOffByDefault()
        {
            Assert.IsFalse(new FlashLfqParameters().WriteVerbosePeaks);

            FlashLfqResults results = Quantify(writeVerbosePeaks: false, out SpectraFileInfo mzml);
            Assert.IsFalse(results.WriteVerbosePeaks);

            ChromatographicPeak peak = results.Peaks[mzml].Single();
            Assert.IsTrue(peak.IsotopicEnvelopes.Any());
            Assert.IsFalse(peak.IsotopicEnvelopes.Any(e => e is VerboseIsotopicEnvelope));

            // a peak with no verbose envelopes still produces the two empty columns rather than a broken row
            Assert.AreEqual("\t\t", peak.GetIsotopeInformation());
            Assert.AreEqual(peak.ToString(), peak.ToString(false));
        }

        /// <summary>
        /// With the option on, every envelope keeps the isotope peaks it was summed from, and the summed intensity
        /// the ordinary columns report is unchanged - the option adds output, it does not change quantification.
        /// </summary>
        [Test]
        public static void TestVerbosePeaksRetainIsotopePeaksWithoutChangingQuantification()
        {
            FlashLfqResults plain = Quantify(writeVerbosePeaks: false, out SpectraFileInfo plainMzml);
            FlashLfqResults verbose = Quantify(writeVerbosePeaks: true, out SpectraFileInfo verboseMzml);

            ChromatographicPeak plainPeak = plain.Peaks[plainMzml].Single();
            ChromatographicPeak verbosePeak = verbose.Peaks[verboseMzml].Single();

            // same quantification either way
            Assert.AreEqual(plainPeak.Intensity, verbosePeak.Intensity, 1e-6);
            Assert.AreEqual(plainPeak.IsotopicEnvelopes.Count, verbosePeak.IsotopicEnvelopes.Count);
            CollectionAssert.AreEqual(
                plainPeak.IsotopicEnvelopes.Select(e => e.Intensity).ToArray(),
                verbosePeak.IsotopicEnvelopes.Select(e => e.Intensity).ToArray());

            // but the verbose run kept the peaks behind each envelope
            var verboseEnvelopes = verbosePeak.IsotopicEnvelopes.OfType<VerboseIsotopicEnvelope>().ToList();
            Assert.AreEqual(verbosePeak.IsotopicEnvelopes.Count, verboseEnvelopes.Count);
            Assert.IsTrue(verboseEnvelopes.All(e => e.PeakDictionary.Count > 0));

            // the monoisotopic peak is present, the isotope numbers ascend from it, and each retained peak really
            // does sit at the m/z its isotope number implies
            foreach (VerboseIsotopicEnvelope envelope in verboseEnvelopes)
            {
                double monoisotopicMz = 1350.65681.ToMz(envelope.ChargeState);
                double spacing = Constants.C13MinusC12 / envelope.ChargeState;
                var tolerance = new PpmTolerance(20);

                foreach (KeyValuePair<int, IIndexedPeak> isotopePeak in envelope.PeakDictionary)
                {
                    Assert.IsTrue(tolerance.Within(isotopePeak.Value.M, monoisotopicMz + isotopePeak.Key * spacing),
                        "isotope " + isotopePeak.Key + " sat at " + isotopePeak.Value.M);
                }

                // every retained peak was observed in the scan the envelope belongs to
                Assert.IsTrue(envelope.PeakDictionary.Values
                    .All(p => p.ZeroBasedScanIndex == envelope.IndexedPeak.ZeroBasedScanIndex));
                Assert.AreEqual(envelope.IndexedPeak.RetentionTime, envelope.RetentionTime);
            }
        }

        /// <summary>
        /// The header gains three columns, so every row has to gain exactly three fields. An earlier draft of this
        /// feature emitted a trailing tab from the isotope fragment, which put a fourth empty column on every
        /// verbose row and shifted nothing visibly - the file simply had one more column than its header.
        /// </summary>
        [Test]
        public static void TestVerbosePeaksFileHasAsManyColumnsAsItsHeader()
        {
            FlashLfqResults results = Quantify(writeVerbosePeaks: true, out SpectraFileInfo mzml);
            Assert.IsTrue(results.WriteVerbosePeaks);

            string peaksPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "verbosePeaks.tsv");
            results.WriteResults(peaksPath, null, null, null, true);

            string[] lines = File.ReadAllLines(peaksPath);
            Assert.IsTrue(lines.Length > 1);

            int headerColumns = lines[0].Split('\t').Length;
            Assert.AreEqual(ChromatographicPeak.TabSeparatedHeader.Split('\t').Length + 3, headerColumns);
            CollectionAssert.AreEqual(new[] { "Isotope Peak Intensity", "Isotope Peak m/z", "Isotope Peak RTs" },
                lines[0].Split('\t').Skip(headerColumns - 3).ToArray());

            foreach (string line in lines.Skip(1))
            {
                Assert.AreEqual(headerColumns, line.Split('\t').Length, line);
            }

            // and the three columns actually carry the envelope: one bracketed list per isotope of per charge state,
            // each holding one entry per retention time in the third column
            string[] fields = lines[1].Split('\t');
            string intensities = fields[headerColumns - 3];
            string mzValues = fields[headerColumns - 2];
            string retentionTimes = fields[headerColumns - 1];

            Assert.IsTrue(intensities.Contains("i0z2"), intensities);
            Assert.IsTrue(mzValues.Contains("i0z2"), mzValues);
            Assert.AreEqual(intensities.Count(c => c == '['), mzValues.Count(c => c == '['));

            int retentionTimeCount = retentionTimes.Trim('[', ']', '"').Split(',').Length;
            string firstIsotopeList = intensities.Split(';').First();
            Assert.AreEqual(retentionTimeCount, firstIsotopeList.Split(',').Length);

            File.Delete(peaksPath);
        }

        /// <summary>
        /// The isotope number is computed from the m/z rather than searched for, so a peak that matches no isotope
        /// position is dropped instead of widening the tolerance until something fits. An earlier draft recursed
        /// with a 5 ppm wider tolerance on every miss, which for a peak sitting between two isotope positions
        /// widened without bound.
        /// </summary>
        [Test]
        public static void TestBuildPeakDictionaryDropsPeaksThatMatchNoIsotope()
        {
            const int chargeState = 2;
            const double monoisotopicMass = 1350.65681;
            double monoisotopicMz = monoisotopicMass.ToMz(chargeState);
            double spacing = Constants.C13MinusC12 / chargeState;
            var tolerance = new PpmTolerance(5);

            var peaks = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(monoisotopicMz, 100, 0, 1.0),
                new IndexedMassSpectralPeak(monoisotopicMz + spacing, 50, 0, 1.0),
                new IndexedMassSpectralPeak(monoisotopicMz + 2 * spacing, 25, 0, 1.0),

                // squarely between the second and third isotope positions, so it belongs to neither
                new IndexedMassSpectralPeak(monoisotopicMz + 1.5 * spacing, 10, 0, 1.0)
            };

            var dictionary = VerboseIsotopicEnvelope.BuildPeakDictionary(peaks, monoisotopicMass, chargeState,
                tolerance);

            CollectionAssert.AreEqual(new[] { 0, 1, 2 }, dictionary.Keys.OrderBy(k => k).ToArray());
            Assert.AreEqual(100, dictionary[0].Intensity, 1e-6);
            Assert.AreEqual(50, dictionary[1].Intensity, 1e-6);
            Assert.AreEqual(25, dictionary[2].Intensity, 1e-6);

            // degenerate inputs come back empty rather than throwing or looping
            Assert.AreEqual(0, VerboseIsotopicEnvelope.BuildPeakDictionary(null, monoisotopicMass, chargeState,
                tolerance).Count);
            Assert.AreEqual(0, VerboseIsotopicEnvelope.BuildPeakDictionary(new List<IIndexedPeak>(), monoisotopicMass,
                chargeState, tolerance).Count);
            Assert.AreEqual(0, VerboseIsotopicEnvelope.BuildPeakDictionary(peaks, monoisotopicMass, 0, tolerance)
                .Count);
        }

        /// <summary>
        /// Two peaks can land on the same isotope number only just inside the tolerance. The more intense one wins,
        /// since that is the one peak finding traced.
        /// </summary>
        [Test]
        public static void TestBuildPeakDictionaryKeepsTheMoreIntenseOfTwoPeaksOnOneIsotope()
        {
            const int chargeState = 2;
            const double monoisotopicMass = 1350.65681;
            double monoisotopicMz = monoisotopicMass.ToMz(chargeState);

            var peaks = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(monoisotopicMz, 10, 0, 1.0),
                new IndexedMassSpectralPeak(monoisotopicMz + monoisotopicMz * 1e-6, 90, 0, 1.0)
            };

            var dictionary = VerboseIsotopicEnvelope.BuildPeakDictionary(peaks, monoisotopicMass, chargeState,
                new PpmTolerance(5));

            Assert.AreEqual(1, dictionary.Count);
            Assert.AreEqual(90, dictionary[0].Intensity, 1e-6);
        }
    }
}
