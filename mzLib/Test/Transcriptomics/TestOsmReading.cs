using NUnit.Framework;
using Omics.Fragmentation;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Fragmentation.Oligo;
using Chemistry;

namespace Test.Transcriptomics;

public class TestOsmReading
{
    public static string OsmFilePath => Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "OsmFileForTesting.osmtsv");
    public static IEnumerable<TestCaseData> IonTestCases
    {
        get
        {
            // All expected values from OsmFileForTesting.osmtsv, first row
            yield return new TestCaseData(
                "a14-3", "a", 14, 1469.18619, 4293, -3, -0.06710, -15.21, false, 0.0, FragmentationTerminus.FivePrime
            ).SetName("a14-3");
            yield return new TestCaseData(
                "a18-4", "a", 18, 1428.93092, 3755, -4, -0.07248, -12.67, false, 0.0, FragmentationTerminus.FivePrime
            ).SetName("a18-4");
            yield return new TestCaseData(
                "aBaseLoss3-1", "aBaseLoss", 3, 725.07528, 4101, -1, 0.00026, 0.36, false, 0.0, FragmentationTerminus.FivePrime
            ).SetName("aBaseLoss3-1");
            yield return new TestCaseData(
                "dWaterLoss2-1", "dWaterLoss", 2, 611.04351, 26007, -1, 0.00018, 0.29, false, 0.0, FragmentationTerminus.FivePrime
            ).SetName("dWaterLoss2-1");
            yield return new TestCaseData(
                "w2-1", "w", 2, 628.06982, 22475, -1, -0.00005, -0.08, false, 0.0, FragmentationTerminus.ThreePrime
            ).SetName("w2-1");
            yield return new TestCaseData(
                "y2-1", "y", 2, 548.10419, 19195, -1, 0.00064, 1.17, false, 0.0, FragmentationTerminus.ThreePrime
            ).SetName("y2-1");
            yield return new TestCaseData(
                "yWaterLoss3-1", "yWaterLoss", 3, 875.13937, 2385, -1, -0.00104, -1.19, false, 0.0, FragmentationTerminus.ThreePrime
            ).SetName("yWaterLoss3-1");
        }
    }

    [Test]
    public static void LoadsWithoutCrashing_OsmSpecific()
    {
        List<string> errors = [];
        List<OsmFromTsv> results = [];
        Assert.DoesNotThrow(() => results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors));
        Assert.That(errors.Count, Is.EqualTo(0));
        Assert.That(results.Count, Is.EqualTo(6));
    }

    [Test]
    public static void LoadsWithoutCrashing_Generic()
    {
        List<string> errors = [];
        List<SpectrumMatchFromTsv> results = [];
        Assert.DoesNotThrow(() => results = SpectrumMatchTsvReader.ReadTsv<SpectrumMatchFromTsv>(OsmFilePath, out errors));
        Assert.That(errors.Count, Is.EqualTo(0));
        Assert.That(results.Count, Is.EqualTo(6));
    }

    [Test, TestCaseSource(nameof(IonTestCases))]
    public static void MatchedFragmentIonProperties_AreCorrectlySet_Case(
        string annotation,
        string productType,
        int fragmentNumber,
        double mz,
        double intensity,
        int charge,
        double massErrorDa,
        double massErrorPpm,
        bool isInternalFragment,
        double neutralLoss, 
        FragmentationTerminus terminus
    )
    {
        List<string> errors = [];
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);

        Assert.That(results.Count, Is.GreaterThan(0), "No results loaded from OSM TSV file.");

        var first = results.First();
        var ion = first.MatchedIons.FirstOrDefault(i => i.Annotation.StartsWith(annotation));
        Assert.That(ion, Is.Not.Null, $"Ion '{annotation}' not found.");

        // MatchedFragmentIon properties
        Assert.That(Math.Round(ion.Mz, 5), Is.EqualTo(mz).Within(1e-4), $"Mz not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.Intensity, 0), Is.EqualTo(intensity).Within(1), $"Intensity not set correctly for {annotation}.");
        Assert.That(ion.Charge, Is.EqualTo(charge), $"Charge not set correctly for {annotation}.");
        Assert.That(ion.Annotation.StartsWith(annotation), $"Annotation not set correctly for {annotation}.");
        Assert.That(ion.IsInternalFragment, Is.EqualTo(isInternalFragment), $"IsInternalFragment not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.MassErrorDa, 5), Is.EqualTo(massErrorDa).Within(1e-4), $"MassErrorDa not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.MassErrorPpm, 2), Is.EqualTo(massErrorPpm).Within(0.1), $"MassErrorPpm not set correctly for {annotation}.");

        // NeutralTheoreticalProduct (Product) properties
        Assert.That(ion.NeutralTheoreticalProduct.ProductType.ToString(), Is.EqualTo(productType), $"ProductType not set correctly for {annotation}.");
        Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(fragmentNumber), $"FragmentNumber not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.NeutralTheoreticalProduct.NeutralMass, 5), Is.EqualTo(mz.ToMass(ion.Charge)).Within(0.1), $"NeutralMass not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.NeutralTheoreticalProduct.MonoisotopicMass, 5), Is.EqualTo(mz.ToMass(ion.Charge)).Within(0.1), $"MonoisotopicMass not set correctly for {annotation}.");
        Assert.That(Math.Round(ion.NeutralTheoreticalProduct.NeutralLoss, 2), Is.EqualTo(neutralLoss).Within(0.01), $"NeutralLoss not set correctly for {annotation}.");
        Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(terminus), $"Terminus not set correctly for {annotation}.");
    }  
}
