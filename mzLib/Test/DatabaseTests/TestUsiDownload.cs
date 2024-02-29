using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Newtonsoft.Json;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    internal class TestUsiDownload
    {
        public const string Prefix = "mzspec";
        public const string Identifier = "PXD000561";
        public const string Run = "Adult_Frontalcortex_bRP_Elite_85_f09";
        public const string TypeFlag = "scan";
        public const string Index = "17555";
        public const string Interpretation = "VLHPLEGAVVIIFK/2";

        [Test]
        public static void TestTryGetSpectraWithoutInterpretation()
        {
            string usi = string.Join(':', new List<string> { Prefix, Identifier, Run, TypeFlag, Index });
            Assert.True(UsiLoader.TryGetSpectra(usi, out var spectrumList));
            Assert.NotNull(spectrumList);
            Assert.AreEqual(1, spectrumList.Count);

            UsiSpectrumFromJSON usiSpectrum = spectrumList.First();
            MzSpectrum mzSpectrum = usiSpectrum.GetMzSpectrum();
            Assert.NotNull(mzSpectrum);
            Assert.AreEqual(350.2189, mzSpectrum.XofPeakWithHighestY);

            // Test getter with indivdual usi elements as arguments
            Assert.True(UsiLoader.TryGetSpectra(Identifier, Run, TypeFlag, Index, out spectrumList));
            Assert.NotNull(spectrumList);
            Assert.AreEqual(1, spectrumList.Count);

            usiSpectrum = spectrumList.First();
            mzSpectrum = usiSpectrum.GetMzSpectrum();
            Assert.NotNull(mzSpectrum);
            Assert.AreEqual(350.2189, mzSpectrum.XofPeakWithHighestY);
            Assert.AreEqual(767.9700, usiSpectrum.GetIsolationWindowTargetMz());
        }

        [Test]
        public static void TestTryGetSpectraWithInterpretation()
        {
            Assert.True(UsiLoader.TryGetSpectra(Identifier, Run, TypeFlag, Index, Interpretation, out var spectrumList));
            Assert.NotNull(spectrumList);
            Assert.AreEqual(1, spectrumList.Count);

            UsiSpectrumFromJSON usiSpectrum = spectrumList.First();
            MzSpectrum mzSpectrum = usiSpectrum.GetMzSpectrum();
            // Note: Checking this spectrum in an online viewer lists the b3-ion as the most intense, with a theoretical mass of 350.2187. However, fragment masses are theoretical, calculated from the peptide sequence. Not experimental
            Assert.AreEqual(350.2189, mzSpectrum.XofPeakWithHighestY); 
            Assert.True(usiSpectrum.TryGetBaseSequence(out string sequence));
            Assert.AreEqual("VLHPLEGAVVIIFK", sequence);

        }


        [Test]
        public static void TestUsiDeserializer()
        {
            var usiSpectraList = JsonConvert.DeserializeObject<List<UsiSpectrumFromJSON>>(ExampleResponse);
            var usiSpectrum = usiSpectraList.First();

            Assert.True(usiSpectrum.TryGetBaseSequence(out string sequence));
            Assert.AreEqual("VLHPLEGAVVIIFK", sequence);
            Assert.AreEqual(2, usiSpectrum.GetMzSpectrum().Size);
        }

        public const string ExampleResponse = @"[
  {
    ""attributes"": [
      {
        ""accession"": ""MS:1008025"",
        ""name"": ""scan number"",
        ""value"": ""17555""
      },
      {
        ""accession"": ""MS:1000827"",
        ""name"": ""isolation window target m/z"",
        ""value"": ""767.9700""
      },
      {
        ""accession"": ""MS:1000041"",
        ""name"": ""charge state"",
        ""value"": ""2""
      },
      {
        ""accession"": ""MS:1003061"",
        ""name"": ""spectrum name"",
        ""value"": ""VLHPLEGAVVIIFK/2""
      },
      {
        ""accession"": ""MS:1000888"",
        ""name"": ""unmodified peptide sequence"",
        ""value"": ""VLHPLEGAVVIIFK""
      }
    ],
    ""intensities"": [
      39316.4648,
      319.6931
    ],
    ""mzs"": [
      110.0712,
      111.0682
    ]
  }
]
";
    }
}
