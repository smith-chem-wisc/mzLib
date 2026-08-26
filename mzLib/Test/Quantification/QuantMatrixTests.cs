using MassSpectrometry;
using NUnit.Framework;
using Omics.SpectralMatch;
using Quantification;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Test.Omics;

namespace Test.Quantification;

/// <summary>
/// Round-trip and boundary coverage for <see cref="QuantMatrix{T}"/>, the structure every
/// quantification strategy reads and writes.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class QuantMatrixTests
{
    private static MockSpectralMatch Psm(string file, string sequence, int scan) =>
        new(file, sequence, sequence, 50.0, scan);

    private static (SpectralMatchMatrix matrix, MockSpectralMatch a, MockSpectralMatch b, MockSpectralMatch c)
        BuildMatrix()
    {
        var psm1 = Psm("Sample1.raw", "PEPTIDEK", 100);
        var psm2 = Psm("Sample1.raw", "ANOTHERPEPTIDEK", 200);
        var psm3 = Psm("Sample2.raw", "THIRDPEPTIDEK", 150);

        var samples = new List<ISampleInfo>
        {
            new SpectraFileInfo("Sample1.raw", "Control", biorep: 1, techrep: 1, fraction: 1),
            new SpectraFileInfo("Sample2.raw", "Treatment", biorep: 1, techrep: 1, fraction: 1)
        };

        var matrix = new SpectralMatchMatrix(
            new List<ISpectralMatch> { psm1, psm2, psm3 }, samples);

        return (matrix, psm1, psm2, psm3);
    }

    [Test]
    public void SetRow_ThenGetRow_RoundTripsByKeyAndByIndex()
    {
        var (matrix, psm1, psm2, psm3) = BuildMatrix();

        matrix.SetRow(psm1, new[] { 1000.0, 500.0 });
        matrix.SetRow(psm2, new[] { 2000.0, 1500.0 });
        matrix.SetRow(psm3, new[] { 0.0, 3000.0 });   // psm3 not detected in Sample1

        Assert.Multiple(() =>
        {
            Assert.That(matrix.GetRow(psm1), Is.EqualTo(new[] { 1000.0, 500.0 }));
            Assert.That(matrix.GetRow(psm2), Is.EqualTo(new[] { 2000.0, 1500.0 }));
            Assert.That(matrix.GetRow(psm3), Is.EqualTo(new[] { 0.0, 3000.0 }));
            Assert.That(matrix.GetRow(0), Is.EqualTo(matrix.GetRow(psm1)));
            Assert.That(matrix.RowCount, Is.EqualTo(3));
            Assert.That(matrix.ColumnCount, Is.EqualTo(2));
        });
    }

    [Test]
    public void SetColumn_WritesEveryRowOfThatSample()
    {
        var (matrix, psm1, psm2, psm3) = BuildMatrix();
        var sample2 = matrix.ColumnKeys[1];

        matrix.SetColumn(sample2, new[] { 11.0, 22.0, 33.0 });

        Assert.Multiple(() =>
        {
            Assert.That(matrix.GetRow(psm1)[1], Is.EqualTo(11.0));
            Assert.That(matrix.GetRow(psm2)[1], Is.EqualTo(22.0));
            Assert.That(matrix.GetRow(psm3)[1], Is.EqualTo(33.0));
            Assert.That(matrix.GetRow(psm1)[0], Is.EqualTo(0.0), "the other column is untouched");
        });
    }

    /// <summary>
    /// Strategies rent buffers from an ArrayPool, which returns arrays at least as long as requested.
    /// A longer-than-needed buffer must be accepted, and only the first ColumnCount entries read.
    /// </summary>
    [Test]
    public void SetRow_AcceptsAnOversizedBuffer_AndReadsOnlyTheColumnsItNeeds()
    {
        var (matrix, psm1, _, _) = BuildMatrix();

        matrix.SetRow(psm1, new[] { 1.0, 2.0, 999.0, 998.0 });

        Assert.That(matrix.GetRow(psm1), Is.EqualTo(new[] { 1.0, 2.0 }));
    }

    [Test]
    public void SetRow_WithUnknownKey_Throws()
    {
        var (matrix, _, _, _) = BuildMatrix();
        var stranger = Psm("Sample3.raw", "UNSEENK", 1);

        Assert.That(() => matrix.SetRow(stranger, new[] { 1.0, 2.0 }),
            Throws.ArgumentException.With.Message.Contains("Row key not found"));
    }

    [Test]
    public void SetRow_WithTooFewValues_Throws()
    {
        var (matrix, psm1, _, _) = BuildMatrix();

        Assert.That(() => matrix.SetRow(psm1, new[] { 1.0 }),
            Throws.ArgumentException.With.Message.Contains("smaller than number of columns"));
    }

    [Test]
    public void GetRow_WithUnknownKeyOrOutOfRangeIndex_Throws()
    {
        var (matrix, _, _, _) = BuildMatrix();
        var stranger = Psm("Sample3.raw", "UNSEENK", 1);

        Assert.Multiple(() =>
        {
            Assert.That(() => matrix.GetRow(stranger), Throws.ArgumentException);
            Assert.That(() => matrix.GetRow(-1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => matrix.GetRow(matrix.RowCount), Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void SetColumn_WithUnknownKeyOrTooFewValues_Throws()
    {
        var (matrix, _, _, _) = BuildMatrix();
        var stranger = new SpectraFileInfo("Sample9.raw", "Other", 1, 1, 1);

        Assert.Multiple(() =>
        {
            Assert.That(() => matrix.SetColumn(stranger, new[] { 1.0, 2.0, 3.0 }),
                Throws.ArgumentException.With.Message.Contains("Column key not found"));
            Assert.That(() => matrix.SetColumn(matrix.ColumnKeys[0], new[] { 1.0 }),
                Throws.ArgumentException.With.Message.Contains("smaller than number of rows"));
        });
    }

    /// <summary>A freshly constructed matrix reads as all zeros rather than uninitialized memory.</summary>
    [Test]
    public void NewMatrix_IsZeroFilled()
    {
        var (matrix, psm1, _, _) = BuildMatrix();
        Assert.That(matrix.GetRow(psm1), Is.EqualTo(new[] { 0.0, 0.0 }));
    }
}
