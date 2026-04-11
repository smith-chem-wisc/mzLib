using System;
using System.IO;
using NUnit.Framework;
using TopDownSimulator.Extraction;

namespace Test.TopDownSimulator;

[TestFixture]
public class MmResultLoaderTests
{
    [Test]
    public void LoadsMetaMorpheusRecordsFromPsmTsv()
    {
        string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"mm-loader-{Guid.NewGuid():N}.psmtsv");
        try
        {
            File.WriteAllText(path,
                "File Name\tScan Number\tPrecursor Scan Number\tPrecursor Charge\tPrecursor MZ\tPrecursor Mass\tBase Sequence\tFull Sequence\tMonoisotopic Mass\tScore\tDecoy/Contaminant/Target\tQValue\tScan Retention Time\tAccession\n" +
                "sample.raw\t100\t95\t12\t1000\t11988\tPEPTIDE\tPEPTIDE\t12345.67\t42\tT\t0.01\t25.5\tP12345\n");

            var loader = new MmResultLoader();
            var results = loader.Load(path);

            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].FileNameWithoutExtension, Is.EqualTo("sample"));
            Assert.That(results[0].Ms2ScanNumber, Is.EqualTo(100));
            Assert.That(results[0].PrecursorScanNumber, Is.EqualTo(95));
            Assert.That(results[0].PrecursorCharge, Is.EqualTo(12));
            Assert.That(results[0].MonoisotopicMass, Is.EqualTo(12345.67).Within(1e-9));
            Assert.That(results[0].RetentionTime, Is.EqualTo(25.5).Within(1e-9));
            Assert.That(results[0].Identifier, Does.Contain("P12345"));
        }
        finally
        {
            if (File.Exists(path))
                File.Delete(path);
        }
    }
}
