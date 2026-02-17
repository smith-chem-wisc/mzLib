using NUnit.Framework;
using Omics.Fragmentation;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Fragmentation.Oligo;
using Chemistry;
using Transcriptomics;
using Transcriptomics.Digestion;

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
        List<string> errors = new();
        List<OsmFromTsv> results = new();
        Assert.DoesNotThrow(() => results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors));
        Assert.That(errors.Count, Is.EqualTo(0));
        Assert.That(results.Count, Is.EqualTo(6));
    }

    [Test]
    public static void LoadsWithoutCrashing_Generic()
    {
        List<string> errors = new();
        List<SpectrumMatchFromTsv> results = new();
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
        List<string> errors = new();
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

    [Test]
    public static void TerminusProperties_AssignedCorrectly_WhenPreviousAndNextResiduesAreTerminal()
    {
        // Test that FivePrimeTerminus and ThreePrimeTerminus are correctly assigned based on PreviousResidue/NextResidue
        List<string> errors = new();
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);

        Assert.That(results.Count, Is.GreaterThan(0), "No results loaded from OSM TSV file.");

        foreach (var result in results)
        {
            // FivePrimeTerminus: if PreviousResidue == "-" => NucleicAcid.DefaultFivePrimeTerminus, else Rnase.DefaultFivePrimeTerminus
            string actualFive = result.FivePrimeTerminus.ThisChemicalFormula.Formula;
            string expectedFive = result.PreviousResidue == "-"
                ? NucleicAcid.DefaultFivePrimeTerminus.Formula
                : Rnase.DefaultFivePrimeTerminus.ThisChemicalFormula.Formula;
            Assert.That(actualFive, Is.EqualTo(expectedFive), 
                $"FivePrimeTerminus mismatch for scan {result.Ms2ScanNumber}. Expected: {expectedFive}, Got: {actualFive}");

            // ThreePrimeTerminus: if NextResidue == "-" => NucleicAcid.DefaultThreePrimeTerminus, else Rnase.DefaultThreePrimeTerminus
            string actualThree = result.ThreePrimeTerminus.ThisChemicalFormula.Formula;
            string expectedThree = result.NextResidue == "-"
                ? NucleicAcid.DefaultThreePrimeTerminus.Formula
                : Rnase.DefaultThreePrimeTerminus.ThisChemicalFormula.Formula;
            Assert.That(actualThree, Is.EqualTo(expectedThree), 
                $"ThreePrimeTerminus mismatch for scan {result.Ms2ScanNumber}. Expected: {expectedThree}, Got: {actualThree}");
        }
    }

    [Test]
    public static void TerminusProperties_NotNull_ForAllResults()
    {
        // Ensure terminus properties are always set and never null
        List<string> errors = new();
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);

        Assert.That(results.Count, Is.GreaterThan(0));

        foreach (var result in results)
        {
            Assert.That(result.FivePrimeTerminus, Is.Not.Null, 
                $"FivePrimeTerminus is null for scan {result.Ms2ScanNumber}");
            Assert.That(result.ThreePrimeTerminus, Is.Not.Null, 
                $"ThreePrimeTerminus is null for scan {result.Ms2ScanNumber}");
            Assert.That(result.FivePrimeTerminus.ThisChemicalFormula, Is.Not.Null,
                $"FivePrimeTerminus.ThisChemicalFormula is null for scan {result.Ms2ScanNumber}");
            Assert.That(result.ThreePrimeTerminus.ThisChemicalFormula, Is.Not.Null,
                $"ThreePrimeTerminus.ThisChemicalFormula is null for scan {result.Ms2ScanNumber}");
        }
    }

    [Test]
    public static void DisambiguatingConstructor_PreservesTerminusProperties()
    {
        // Test the disambiguating constructor (OsmFromTsv(OsmFromTsv osm, string fullSequence, ...))
        List<string> errors = new();
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);
        var originalOsm = results.First();

        // Test 1: Constructor without explicit terminus parameters (should copy from original)
        var clonedOsm = new OsmFromTsv(originalOsm, originalOsm.FullSequence);
        
        Assert.That(clonedOsm.FivePrimeTerminus.ThisChemicalFormula.Formula, 
            Is.EqualTo(originalOsm.FivePrimeTerminus.ThisChemicalFormula.Formula),
            "Cloned OSM should preserve FivePrimeTerminus when not explicitly provided");
        Assert.That(clonedOsm.ThreePrimeTerminus.ThisChemicalFormula.Formula, 
            Is.EqualTo(originalOsm.ThreePrimeTerminus.ThisChemicalFormula.Formula),
            "Cloned OSM should preserve ThreePrimeTerminus when not explicitly provided");
    }

    [Test]
    public static void DisambiguatingConstructor_OverridesTerminusProperties_WhenProvided()
    {
        // Test the disambiguating constructor with explicit terminus parameters
        List<string> errors = new();
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);
        var originalOsm = results.First();

        // Create custom terminus formulas
        var customFivePrime = ChemicalFormula.ParseFormula("H2O");
        var customThreePrime = ChemicalFormula.ParseFormula("PO4");

        // Test 2: Constructor with explicit terminus parameters (should override original)
        var customOsm = new OsmFromTsv(originalOsm, originalOsm.FullSequence, 
            fivePrimeTerm: customFivePrime, threePrimeTerm: customThreePrime);
        
        Assert.That(customOsm.FivePrimeTerminus.ThisChemicalFormula.Formula, 
            Is.EqualTo(customFivePrime.Formula),
            "Custom FivePrimeTerminus should override original");
        Assert.That(customOsm.ThreePrimeTerminus.ThisChemicalFormula.Formula, 
            Is.EqualTo(customThreePrime.Formula),
            "Custom ThreePrimeTerminus should override original");
    }

    [Test]
    public static void OsmFromTsvFile_LoadsCorrectly()
    {
        // Test loading via OsmFromTsvFile
        var osmFile = new OsmFromTsvFile(OsmFilePath);
        osmFile.LoadResults();

        Assert.That(osmFile.Results, Is.Not.Null);
        Assert.That(osmFile.Results.Count, Is.EqualTo(6));
        Assert.That(osmFile.FileType, Is.EqualTo(SupportedFileType.osmtsv));

        // Verify terminus properties are set for all results
        foreach (var result in osmFile.Results)
        {
            Assert.That(result.FivePrimeTerminus, Is.Not.Null);
            Assert.That(result.ThreePrimeTerminus, Is.Not.Null);
        }
    }

    [Test]
    public static void TerminalCapFormulasAreConsistent()
    {
        string expectedNucleicAcidFivePrime = "O-3P-1";
        Assert.That(NucleicAcid.DefaultFivePrimeTerminus.Formula, Is.EqualTo(expectedNucleicAcidFivePrime),
            "NucleicAcid.DefaultFivePrimeTerminus formula has changed unexpectedly.");

        string expectedNucleicAcidThreePrime = "HO";
        Assert.That(NucleicAcid.DefaultThreePrimeTerminus.Formula, Is.EqualTo(expectedNucleicAcidThreePrime),
            "NucleicAcid.DefaultThreePrimeTerminus formula has changed unexpectedly.");

        string expectedRnaseFivePrime = "O-3P-1";
        Assert.That(Rnase.DefaultFivePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(expectedRnaseFivePrime),
            "Rnase.DefaultFivePrimeTerminus formula has changed unexpectedly.");

        string expectedRnaseThreePrime = "H2O4P";
        Assert.That(Rnase.DefaultThreePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(expectedRnaseThreePrime),
            "Rnase.DefaultThreePrimeTerminus formula has changed unexpectedly.");
    }

    [Test]
    public static void TerminusProperties_MatchExpectedFormulas()
    {
        // Verify the specific chemical formulas match expected values
        List<string> errors = new();
        var results = SpectrumMatchTsvReader.ReadOsmTsv(OsmFilePath, out errors);

        // Expected formulas from the static properties
        IHasChemicalFormula expectedNucleicAcidFivePrime = NucleicAcid.DefaultFivePrimeTerminus;
        IHasChemicalFormula expectedNucleicAcidThreePrime = NucleicAcid.DefaultThreePrimeTerminus;
        IHasChemicalFormula expectedRnaseFivePrime = Rnase.DefaultFivePrimeTerminus;
        IHasChemicalFormula expectedRnaseThreePrime = Rnase.DefaultThreePrimeTerminus;

        foreach (var result in results)
        {
            if (result.PreviousResidue == "-")
            {
                Assert.That(result.FivePrimeTerminus.ThisChemicalFormula, 
                    Is.EqualTo(expectedNucleicAcidFivePrime),
                    $"Terminal oligo (PreviousResidue='-') should have NucleicAcid.DefaultFivePrimeTerminus for scan {result.Ms2ScanNumber}");
            }
            else
            {
                Assert.That(result.FivePrimeTerminus.ThisChemicalFormula, 
                    Is.EqualTo(expectedRnaseFivePrime),
                    $"Internal oligo (PreviousResidue!='-') should have Rnase.DefaultFivePrimeTerminus for scan {result.Ms2ScanNumber}");
            }

            if (result.NextResidue == "-")
            {
                Assert.That(result.ThreePrimeTerminus.ThisChemicalFormula, 
                    Is.EqualTo(expectedNucleicAcidThreePrime),
                    $"Terminal oligo (NextResidue='-') should have NucleicAcid.DefaultThreePrimeTerminus for scan {result.Ms2ScanNumber}");
            }
            else
            {
                Assert.That(result.ThreePrimeTerminus.ThisChemicalFormula, 
                    Is.EqualTo(expectedRnaseThreePrime),
                    $"Internal oligo (NextResidue!='-') should have Rnase.DefaultThreePrimeTerminus for scan {result.Ms2ScanNumber}");
            }
        }
    }

    [Test]
    public static void SpectrumMatchTsvReader_ParsesTerminusHeadersCorrectly()
    {
        // Verify that SpectrumMatchTsvReader.ParseHeader correctly identifies terminus columns
        string testHeader = "File Name\tScan Number\tFull Sequence\t5'-Terminus\t3'-Terminus\tPrevious Residue\tNext Residue";
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(testHeader);

        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.FivePrimeTerminus), Is.True,
            "ParseHeader should recognize 5'-Terminus column");
        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.ThreePrimeTerminus), Is.True,
            "ParseHeader should recognize 3'-Terminus column");
        
        Assert.That(parsedHeader[SpectrumMatchFromTsvHeader.FivePrimeTerminus], Is.EqualTo(3),
            "5'-Terminus should be at index 3");
        Assert.That(parsedHeader[SpectrumMatchFromTsvHeader.ThreePrimeTerminus], Is.EqualTo(4),
            "3'-Terminus should be at index 4");
    }

    [Test]
    public static void SpectrumMatchTsvReader_HandlesAbsentTerminusHeaders()
    {
        // Verify that ParseHeader handles missing terminus columns correctly
        string testHeader = "File Name\tScan Number\tFull Sequence\tPrevious Residue\tNext Residue";
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(testHeader);

        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.FivePrimeTerminus), Is.True,
            "ParseHeader should include FivePrimeTerminus key even when column absent");
        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.ThreePrimeTerminus), Is.True,
            "ParseHeader should include ThreePrimeTerminus key even when column absent");
        
        Assert.That(parsedHeader[SpectrumMatchFromTsvHeader.FivePrimeTerminus], Is.EqualTo(-1),
            "Missing 5'-Terminus column should have index -1");
        Assert.That(parsedHeader[SpectrumMatchFromTsvHeader.ThreePrimeTerminus], Is.EqualTo(-1),
            "Missing 3'-Terminus column should have index -1");
    }

    [Test]
    public static void OSM_ParsesTerminusCorrectly_NotInOsmFile_FullTranscript()
    {
        var osmLines = File.ReadAllLines(OsmFilePath);
        var header = osmLines[0];
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(header);
        var firstDataLine = osmLines[1];
        var spl = firstDataLine.Split('\t').Select(p => p.Trim('\"')).ToArray();

        // Ensure prev and next residues match a full transcript
        spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]] = "-";
        spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]] = "-";
        firstDataLine = string.Join("\t", spl);

        var osm = new OsmFromTsv(firstDataLine, ['\t'], parsedHeader);
        Assert.That(osm.ThreePrimeTerminus, Is.EqualTo(NucleicAcid.DefaultThreePrimeTerminus), 
            "Full transcript should have NucleicAcid.DefaultThreePrimeTerminus");
        Assert.That(osm.FivePrimeTerminus, Is.EqualTo(NucleicAcid.DefaultFivePrimeTerminus),
            "Full transcript should have NucleicAcid.DefaultFivePrimeTerminus");
    }

    [Test]
    public static void OSM_ParsesTerminusCorrectly_NotInOsmFile_TerminalOligo()
    {
        var osmLines = File.ReadAllLines(OsmFilePath);
        var header = osmLines[0];
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(header);
        var firstDataLine = osmLines[1];
        var spl = firstDataLine.Split('\t').Select(p => p.Trim('\"')).ToArray();

        // Ensure prev and next residues match a 5'-terminal oligo
        spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]] = "-";
        spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]] = "A";
        firstDataLine = string.Join("\t", spl);

        var osm = new OsmFromTsv(firstDataLine, ['\t'], parsedHeader);
        Assert.That(osm.FivePrimeTerminus, Is.EqualTo(NucleicAcid.DefaultFivePrimeTerminus),
            "5'-terminal oligo should have NucleicAcid.DefaultFivePrimeTerminus");
        Assert.That(osm.ThreePrimeTerminus, Is.EqualTo(Rnase.DefaultThreePrimeTerminus),
            "5'-terminal oligo should have Rnase.DefaultThreePrimeTerminus");

        // ensure prev and next residues match a 3'-terminal oligo
        spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]] = "U";
        spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]] = "-";
        firstDataLine = string.Join("\t", spl);

        osm = new OsmFromTsv(firstDataLine, ['\t'], parsedHeader);
        Assert.That(osm.FivePrimeTerminus, Is.EqualTo(Rnase.DefaultFivePrimeTerminus),
            "5'-terminal oligo should have Rnase.DefaultFivePrimeTerminus");
        Assert.That(osm.ThreePrimeTerminus, Is.EqualTo(NucleicAcid.DefaultThreePrimeTerminus),
            "3'-terminal oligo should have NucleicAcid.DefaultThreePrimeTerminus");
    }

    [Test]
    public static void OSM_ParsesTerminusCorrectly_NotInOsmFile_InternalOligo()
    {
        var osmLines = File.ReadAllLines(OsmFilePath);
        var header = osmLines[0];
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(header);
        var firstDataLine = osmLines[1];
        var spl = firstDataLine.Split('\t').Select(p => p.Trim('\"')).ToArray();

        // Ensure prev and next residues match a internal oligo
        spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]] = "U";
        spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]] = "A";
        firstDataLine = string.Join("\t", spl);

        var osm = new OsmFromTsv(firstDataLine, ['\t'], parsedHeader);
        Assert.That(osm.ThreePrimeTerminus, Is.EqualTo(Rnase.DefaultThreePrimeTerminus),
            "Internal oligo should have Rnase.DefaultThreePrimeTerminus");
        Assert.That(osm.FivePrimeTerminus, Is.EqualTo(Rnase.DefaultFivePrimeTerminus),
            "Internal oligo should have Rnase.DefaultFivePrimeTerminus");
    }

    [Test]
    public static void OSM_ParsesTerminusCorrectly_InOsmFile_CustomFormulas()
    {
        var osmLines = File.ReadAllLines(OsmFilePath);
        var header = osmLines[0] + "\t5'-Terminus\t3'-Terminus";
        var parsedHeader = SpectrumMatchTsvReader.ParseHeader(header);

        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.FivePrimeTerminus), Is.True,
            "ParseHeader should recognize 5'-Terminus column");
        Assert.That(parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.ThreePrimeTerminus), Is.True,
            "ParseHeader should recognize 3'-Terminus column");

        var dummyFivePrime = ChemicalFormula.ParseFormula("Cr2N7O2");
        var dummyThreePrime = ChemicalFormula.ParseFormula("Y3Ir2Dy6");

        var firstDataLine = osmLines[1] + $"\t\"{dummyFivePrime.Formula}\"\t\"{dummyThreePrime.Formula}\"";

        var osm = new OsmFromTsv(firstDataLine, ['\t'], parsedHeader);
        Assert.That(osm.FivePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(dummyFivePrime.Formula),
            "Custom 5'-Terminus formula not parsed correctly.");
        Assert.That(osm.ThreePrimeTerminus.ThisChemicalFormula.Formula, Is.EqualTo(dummyThreePrime.Formula),
            "Custom 3'-Terminus formula not parsed correctly.");
    }
}
