using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;

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
        public static void TestTryGetSpectra()
        {
            string usi = string.Join(':', new List<string> { Prefix, Identifier, Run, TypeFlag, Index });
            Assert.True(UsiLoader.TryGetSpectrum(usi, out var spectrum));
            Assert.NotNull(spectrum);
        }

        [Test]
        public static void TestTryGetSpectraWithResults()
        {
            string usi = string.Join(':', new List<string> { Prefix, Identifier, Run, TypeFlag, Index, Interpretation });
            Assert.True(UsiLoader.TryGetSpectrum(usi, out var spectrum, out var results));
            Assert.NotNull(spectrum);
            Assert.NotNull(results);
        }
    }
}
