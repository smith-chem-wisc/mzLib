using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test.Transcriptomics;

public class TestOsmReading
{
    public static string OsmFilePath => Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "OsmFileForTesting.osmtsv");

    [Test]
    public static void LoadsWithoutCrashing()
    {
        List<string> errors = [];
        List<SpectrumMatchFromTsv> results = [];
        Assert.DoesNotThrow(() => SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors));
        Assert.That(errors.Count, Is.EqualTo(0));
    }
}
