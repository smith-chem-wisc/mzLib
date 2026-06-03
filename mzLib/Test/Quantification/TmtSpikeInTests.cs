using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Strategies;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Test.Omics;
using Omics.SpectralMatch;

namespace Test.Quantification;

/// <summary>
/// Unit tests for TMT quantification pipeline components using synthetic data.
/// Tests that require real spike-in data files are in the Development project
/// (Development/QuantificationDevelopment/TmtSpikeInDevelopmentTests.cs).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TmtSpikeInTests
{
    #region Private Helper Classes

    private class TmtTestExperimentalDesign : IExperimentalDesign
    {
        public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

        public TmtTestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
        {
            FileNameSampleInfoDictionary = dict;
        }
    }

    #endregion

    #region Helper Methods

    /// <summary>
    /// Creates a TMT experimental design where all files share the same channel-to-condition mapping.
    /// </summary>
    private static IExperimentalDesign CreateTmtExperimentalDesign(
        string[] fileNames,
        string[] channelLabels,
        string[] conditions,
        int[] bioReps,
        int[] techReps)
    {
        var dict = new Dictionary<string, ISampleInfo[]>();
        for (int f = 0; f < fileNames.Length; f++)
        {
            var samples = new ISampleInfo[channelLabels.Length];
            for (int c = 0; c < channelLabels.Length; c++)
            {
                bool isRef = conditions[c] == "Reference";
                samples[c] = new IsobaricQuantSampleInfo(
                    fileNames[f], conditions[c], bioReps[c], techReps[f], 0, 0,
                    channelLabels[c], 126.0 + c * 0.1, isRef);
            }
            dict[fileNames[f]] = samples;
        }
        return new TmtTestExperimentalDesign(dict);
    }

    #endregion

    /// <summary>
    /// Tests that PivotByFile creates one matrix per file, with correct dimensions and values.
    /// </summary>
    [Test]
    public void PivotByFile_CreatesCorrectPerFileMatrices()
    {
        // Setup: 2 files, 3 channels each
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        var file1Psm1 = new MockSpectralMatch(file1, "PEPTIDEK", "PEPTIDEK", 100.0, 1)
            { Intensities = new double[] { 100, 200, 300 } };
        var file1Psm2 = new MockSpectralMatch(file1, "SEQUENCEK", "SEQUENCEK", 95.0, 2)
            { Intensities = new double[] { 400, 500, 600 } };
        var file2Psm1 = new MockSpectralMatch(file2, "PEPTIDEK", "PEPTIDEK", 90.0, 3)
            { Intensities = new double[] { 700, 800, 900 } };
        var file2Psm2 = new MockSpectralMatch(file2, "SEQUENCEK", "SEQUENCEK", 85.0, 4)
            { Intensities = new double[] { 1000, 1100, 1200 } };

        var allPsms = new List<ISpectralMatch> { file1Psm1, file1Psm2, file2Psm1, file2Psm2 };

        var perFileMatrices = QuantificationEngine.PivotByFile(allPsms, design);

        Assert.That(perFileMatrices.Count, Is.EqualTo(2));
        Assert.That(perFileMatrices.ContainsKey(file1), Is.True);
        Assert.That(perFileMatrices.ContainsKey(file2), Is.True);

        var m1 = perFileMatrices[file1];
        Assert.That(m1.RowKeys.Count, Is.EqualTo(2));
        Assert.That(m1.ColumnKeys.Count, Is.EqualTo(3));

        var f1Row0 = m1.GetRow(0);
        Assert.That(f1Row0[0], Is.EqualTo(100.0));
        Assert.That(f1Row0[1], Is.EqualTo(200.0));
        Assert.That(f1Row0[2], Is.EqualTo(300.0));

        var f1Row1 = m1.GetRow(1);
        Assert.That(f1Row1[0], Is.EqualTo(400.0));
        Assert.That(f1Row1[1], Is.EqualTo(500.0));
        Assert.That(f1Row1[2], Is.EqualTo(600.0));

        var m2 = perFileMatrices[file2];
        Assert.That(m2.RowKeys.Count, Is.EqualTo(2));
        Assert.That(m2.ColumnKeys.Count, Is.EqualTo(3));

        var f2Row0 = m2.GetRow(0);
        Assert.That(f2Row0[0], Is.EqualTo(700.0));
        Assert.That(f2Row0[1], Is.EqualTo(800.0));
        Assert.That(f2Row0[2], Is.EqualTo(900.0));

        var f2Row1 = m2.GetRow(1);
        Assert.That(f2Row1[0], Is.EqualTo(1000.0));
        Assert.That(f2Row1[1], Is.EqualTo(1100.0));
        Assert.That(f2Row1[2], Is.EqualTo(1200.0));
    }

    [Test]
    public void CombinePeptideMatrices_MergesCorrectly()
    {
        var protA = new Protein("AAAAK", "PA");
        var protB = new Protein("LLLLK", "PB");
        var protC = new Protein("SSSSSK", "PC");

        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 4);
        var pepA = protA.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pepB = protB.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pepC = protC.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        var file1Samples = design.FileNameSampleInfoDictionary[file1];
        var file2Samples = design.FileNameSampleInfoDictionary[file2];

        var matrix1 = new PeptideMatrix(
            new List<IBioPolymerWithSetMods> { pepA, pepB },
            file1Samples, design);
        matrix1.SetRow(pepA, new double[] { 100, 200, 300 });
        matrix1.SetRow(pepB, new double[] { 400, 500, 600 });

        var matrix2 = new PeptideMatrix(
            new List<IBioPolymerWithSetMods> { pepA, pepC },
            file2Samples, design);
        matrix2.SetRow(pepA, new double[] { 700, 800, 900 });
        matrix2.SetRow(pepC, new double[] { 1000, 1100, 1200 });

        var perFilePeptideMatrices = new Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>>
        {
            { file1, matrix1 },
            { file2, matrix2 }
        };

        var combined = QuantificationEngine.CombinePeptideMatrices(perFilePeptideMatrices, design);

        Assert.That(combined.RowKeys.Count, Is.EqualTo(3));
        Assert.That(combined.ColumnKeys.Count, Is.EqualTo(6));

        var pepARow = combined.GetRow(pepA);
        Assert.That(pepARow[0], Is.EqualTo(100.0));
        Assert.That(pepARow[1], Is.EqualTo(200.0));
        Assert.That(pepARow[2], Is.EqualTo(300.0));
        Assert.That(pepARow[3], Is.EqualTo(700.0));
        Assert.That(pepARow[4], Is.EqualTo(800.0));
        Assert.That(pepARow[5], Is.EqualTo(900.0));

        var pepBRow = combined.GetRow(pepB);
        Assert.That(pepBRow[0], Is.EqualTo(400.0));
        Assert.That(pepBRow[1], Is.EqualTo(500.0));
        Assert.That(pepBRow[2], Is.EqualTo(600.0));
        Assert.That(pepBRow[3], Is.EqualTo(0.0));
        Assert.That(pepBRow[4], Is.EqualTo(0.0));
        Assert.That(pepBRow[5], Is.EqualTo(0.0));

        var pepCRow = combined.GetRow(pepC);
        Assert.That(pepCRow[0], Is.EqualTo(0.0));
        Assert.That(pepCRow[1], Is.EqualTo(0.0));
        Assert.That(pepCRow[2], Is.EqualTo(0.0));
        Assert.That(pepCRow[3], Is.EqualTo(1000.0));
        Assert.That(pepCRow[4], Is.EqualTo(1100.0));
        Assert.That(pepCRow[5], Is.EqualTo(1200.0));
    }

    [Test]
    public void RunTmt_FullPipeline_ProducesCorrectProteinMatrix()
    {
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2 }, channels, conditions, bioReps, techReps);

        var protein1 = new Protein("PEPTIDEK", "P001");
        var protein2 = new Protein("SEQUENCEK", "P002");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();
        var pep2 = protein2.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1, pep2 };

        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var group2 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein2 },
            new HashSet<IBioPolymerWithSetMods> { pep2 },
            new HashSet<IBioPolymerWithSetMods> { pep2 });
        var proteinGroups = new List<IBioPolymerGroup> { group1, group2 };

        var sm1 = new MockSpectralMatch(file1, pep1.FullSequence, pep1.BaseSequence, 100.0, 1, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };
        var sm2 = new MockSpectralMatch(file1, pep2.FullSequence, pep2.BaseSequence, 100.0, 2, new[] { pep2 })
            { Intensities = new double[] { 400, 500, 600 } };
        var sm3 = new MockSpectralMatch(file2, pep1.FullSequence, pep1.BaseSequence, 100.0, 3, new[] { pep1 })
            { Intensities = new double[] { 1000, 2000, 3000 } };
        var sm4 = new MockSpectralMatch(file2, pep2.FullSequence, pep2.BaseSequence, 100.0, 4, new[] { pep2 })
            { Intensities = new double[] { 4000, 5000, 6000 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3, sm4 };

        var parameters = QuantificationParameters.GetSimpleParameters();
        parameters.WriteRawInformation = false;
        parameters.WritePeptideInformation = false;
        parameters.WriteProteinInformation = false;

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);
        engine.Run(out var proteinMatrix);

        Assert.That(proteinMatrix.RowKeys.Count, Is.EqualTo(2));
        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(6));

        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(100.0));
        Assert.That(p1Row[1], Is.EqualTo(200.0));
        Assert.That(p1Row[2], Is.EqualTo(300.0));
        Assert.That(p1Row[3], Is.EqualTo(1000.0));
        Assert.That(p1Row[4], Is.EqualTo(2000.0));
        Assert.That(p1Row[5], Is.EqualTo(3000.0));

        var p2Row = proteinMatrix.GetRow(group2);
        Assert.That(p2Row[0], Is.EqualTo(400.0));
        Assert.That(p2Row[1], Is.EqualTo(500.0));
        Assert.That(p2Row[2], Is.EqualTo(600.0));
        Assert.That(p2Row[3], Is.EqualTo(4000.0));
        Assert.That(p2Row[4], Is.EqualTo(5000.0));
        Assert.That(p2Row[5], Is.EqualTo(6000.0));
    }

    [Test]
    public void RunTmt_WithSumCollapse_CombinesTechnicalReplicates()
    {
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string file3 = "file3.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2, 3 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2, file3 }, channels, conditions, bioReps, techReps);

        var protein1 = new Protein("PEPTIDEK", "P001");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1 };
        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var proteinGroups = new List<IBioPolymerGroup> { group1 };

        var sm1 = new MockSpectralMatch(file1, pep1.FullSequence, pep1.BaseSequence, 100.0, 1, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };
        var sm2 = new MockSpectralMatch(file2, pep1.FullSequence, pep1.BaseSequence, 100.0, 2, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };
        var sm3 = new MockSpectralMatch(file3, pep1.FullSequence, pep1.BaseSequence, 100.0, 3, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3 };

        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new SumCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = string.Empty,
            UseSharedPeptidesForProteinQuant = false,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);
        engine.Run(out var proteinMatrix);

        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(3));

        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(300.0));
        Assert.That(p1Row[1], Is.EqualTo(600.0));
        Assert.That(p1Row[2], Is.EqualTo(900.0));
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Strategy unit tests
    // ──────────────────────────────────────────────────────────────────────────

    [Test]
    public void GlobalMedianNormalization_EqualizesColumnMedians()
    {
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N", "127C" },
            new[] { "CondA", "CondB", "CondC" },
            new[] { 1, 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        var matrix = new QuantMatrix<string>(new[] { "Entity1", "Entity2" }, columns, design);
        matrix.SetRow("Entity1", new double[] { 100, 200, 400 });
        matrix.SetRow("Entity2", new double[] { 100, 200, 400 });

        var norm = new GlobalMedianNormalization();
        var result = norm.NormalizeIntensities(matrix);

        Assert.That(result.Matrix[0, 0], Is.EqualTo(200.0).Within(1e-6), "Col0 row0 should be 200");
        Assert.That(result.Matrix[1, 0], Is.EqualTo(200.0).Within(1e-6), "Col0 row1 should be 200");
        Assert.That(result.Matrix[0, 1], Is.EqualTo(200.0).Within(1e-6), "Col1 row0 should be unchanged");
        Assert.That(result.Matrix[1, 1], Is.EqualTo(200.0).Within(1e-6), "Col1 row1 should be unchanged");
        Assert.That(result.Matrix[0, 2], Is.EqualTo(200.0).Within(1e-6), "Col2 row0 should be 200");
        Assert.That(result.Matrix[1, 2], Is.EqualTo(200.0).Within(1e-6), "Col2 row1 should be 200");

        var matrixWithZero = new QuantMatrix<string>(new[] { "E1", "E2", "E3" }, columns, design);
        matrixWithZero.SetRow("E1", new double[] { 100, 0, 400 });
        matrixWithZero.SetRow("E2", new double[] { 200, 500, 600 });
        matrixWithZero.SetRow("E3", new double[] { 400, 800, 1600 });
        var resultWithZero = norm.NormalizeIntensities(matrixWithZero);
        Assert.That(resultWithZero.Matrix[0, 1], Is.EqualTo(0.0), "Zero should remain zero after normalization");
    }

    [Test]
    public void ReferenceChannelNormalization_ConvertsToRatios()
    {
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N", "131N" },
            new[] { "Reference", "CondA", "Reference" },
            new[] { 1, 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        var matrix = new QuantMatrix<string>(new[] { "Entity1", "Entity2" }, columns, design);
        matrix.SetRow("Entity1", new double[] { 100, 200, 100 });
        matrix.SetRow("Entity2", new double[] { 200, 400, 200 });

        var norm = new ReferenceChannelNormalization();
        var result = norm.NormalizeIntensities(matrix);

        Assert.That(result.Matrix[0, 0], Is.EqualTo(1.0).Within(1e-9), "Entity1 ref ch 126 should be 1.0");
        Assert.That(result.Matrix[0, 2], Is.EqualTo(1.0).Within(1e-9), "Entity1 ref ch 131N should be 1.0");
        Assert.That(result.Matrix[1, 0], Is.EqualTo(1.0).Within(1e-9), "Entity2 ref ch 126 should be 1.0");
        Assert.That(result.Matrix[1, 2], Is.EqualTo(1.0).Within(1e-9), "Entity2 ref ch 131N should be 1.0");

        Assert.That(result.Matrix[0, 1], Is.EqualTo(2.0).Within(1e-9), "Entity1 non-ref ch 127N should be 2.0");
        Assert.That(result.Matrix[1, 1], Is.EqualTo(2.0).Within(1e-9), "Entity2 non-ref ch 127N should be 2.0");
    }

    [Test]
    public void MedianRollUp_ComputesMedianPerColumn()
    {
        var design = CreateTmtExperimentalDesign(
            new[] { "file1.raw" },
            new[] { "126", "127N" },
            new[] { "CondA", "CondB" },
            new[] { 1, 1 }, new[] { 1 });
        var columns = design.FileNameSampleInfoDictionary["file1.raw"];

        var inputMatrix = new QuantMatrix<string>(
            new[] { "PSM1", "PSM2", "PSM3" }, columns, design);
        inputMatrix.SetRow("PSM1", new double[] { 100, 200 });
        inputMatrix.SetRow("PSM2", new double[] { 300, 400 });
        inputMatrix.SetRow("PSM3", new double[] { 500, 100 });

        var map = new Dictionary<string, List<int>>
        {
            { "Protein1", new List<int> { 0, 1, 2 } }
        };

        var rollUp = new MedianRollUp();
        var result = rollUp.RollUp(inputMatrix, map);

        Assert.That(result.Matrix[0, 0], Is.EqualTo(300.0).Within(1e-6), "Median of [100,300,500] should be 300");
        Assert.That(result.Matrix[0, 1], Is.EqualTo(200.0).Within(1e-6), "Median of [200,400,100] should be 200");

        var inputMatrixWithZeros = new QuantMatrix<string>(
            new[] { "PSM1", "PSM2", "PSM3" }, columns, design);
        inputMatrixWithZeros.SetRow("PSM1", new double[] { 100, 0 });
        inputMatrixWithZeros.SetRow("PSM2", new double[] { 300, 400 });
        inputMatrixWithZeros.SetRow("PSM3", new double[] { 0, 100 });

        var resultWithZeros = rollUp.RollUp(inputMatrixWithZeros, map);

        Assert.That(resultWithZeros.Matrix[0, 0], Is.EqualTo(200.0).Within(1e-6), "Median of [100,300] (ignoring zeros) should be 200");
        Assert.That(resultWithZeros.Matrix[0, 1], Is.EqualTo(250.0).Within(1e-6), "Median of [400,100] (ignoring zeros) should be 250");
    }

    [Test]
    public void MeanCollapse_AveragesTechnicalReplicates()
    {
        string file1 = "file1.raw";
        string file2 = "file2.raw";
        string file3 = "file3.raw";
        string[] channels = { "126", "127N", "127C" };
        string[] conditions = { "Reference", "CondA", "CondB" };
        int[] bioReps = { 1, 1, 1 };
        int[] techReps = { 1, 2, 3 };

        var design = CreateTmtExperimentalDesign(
            new[] { file1, file2, file3 }, channels, conditions, bioReps, techReps);

        var protein1 = new Protein("PEPTIDEK", "P001");
        var digParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var pep1 = protein1.Digest(digParams, new List<Modification>(), new List<Modification>()).First();

        var allPeptides = new List<IBioPolymerWithSetMods> { pep1 };
        var group1 = new BioPolymerGroup(
            new HashSet<IBioPolymer> { protein1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 },
            new HashSet<IBioPolymerWithSetMods> { pep1 });
        var proteinGroups = new List<IBioPolymerGroup> { group1 };

        var sm1 = new MockSpectralMatch(file1, pep1.FullSequence, pep1.BaseSequence, 100.0, 1, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };
        var sm2 = new MockSpectralMatch(file2, pep1.FullSequence, pep1.BaseSequence, 100.0, 2, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };
        var sm3 = new MockSpectralMatch(file3, pep1.FullSequence, pep1.BaseSequence, 100.0, 3, new[] { pep1 })
            { Intensities = new double[] { 100, 200, 300 } };

        var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3 };

        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new MeanCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = string.Empty,
            UseSharedPeptidesForProteinQuant = false,
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };

        var engine = new QuantificationEngine(parameters, design, spectralMatches, allPeptides, proteinGroups);
        engine.Run(out var proteinMatrix);

        Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(3),
            "Expected 3 columns after MeanCollapse (Reference_1, CondA_1, CondB_1)");

        var p1Row = proteinMatrix.GetRow(group1);
        Assert.That(p1Row[0], Is.EqualTo(100.0).Within(1e-6), "Reference_1: mean(100,100,100)=100");
        Assert.That(p1Row[1], Is.EqualTo(200.0).Within(1e-6), "CondA_1: mean(200,200,200)=200");
        Assert.That(p1Row[2], Is.EqualTo(300.0).Within(1e-6), "CondB_1: mean(300,300,300)=300");
    }
}
