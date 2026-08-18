using System;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.CZE;
using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.RetentionTimePrediction;

[TestFixture]
public class CZERetentionTimePredictorTests
{
    private sealed class StubRetentionPredictable : IRetentionPredictable
    {
        public string BaseSequence { get; init; } = "PEPTIDE";
        public string FullSequence { get; init; } = "PEPTIDE";
        public string FullSequenceWithMassShifts { get; init; } = "PEPTIDE";
        public double MonoisotopicMass { get; init; } = 1;
    }

    [Test]
    public void Constructor_RejectsNonPositiveColumnLength()
    {
        Assert.Throws<ArgumentException>(() =>
            new CZERetentionTimePredictor(SequenceConversionHandlingMode.UsePrimarySequence, 0, 300000));
    }

    [Test]
    public void Constructor_RejectsNonPositiveVoltage()
    {
        Assert.Throws<ArgumentException>(() =>
            new CZERetentionTimePredictor(SequenceConversionHandlingMode.UsePrimarySequence, 1, 0));
    }

    [Test]
    public void PredictRetentionTime_WithInvalidMass_ReturnsInvalidMassFailure()
    {
        var predictor = new CZERetentionTimePredictor();
        var peptide = new StubRetentionPredictable { MonoisotopicMass = 0 };

        var predicted = predictor.PredictRetentionTimeEquivalent(peptide, out var failureReason);

        Assert.That(predicted, Is.Null);
        Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.InvalidMass));
    }

    [Test]
    public void GetFormattedSequence_ReturnsBaseSequence()
    {
        var predictor = new CZERetentionTimePredictor();
        var peptide = new StubRetentionPredictable { BaseSequence = "ABC" };

        var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);

        Assert.That(formatted, Is.EqualTo("ABC"));
        Assert.That(failureReason, Is.Null);
    }

    [Test]
    public void ElutionAndMobility_InvalidInputs_ReturnMinusOne()
    {
        var predictor = new CZERetentionTimePredictor();

        Assert.That(predictor.TheoreticalElutionTime(0), Is.EqualTo(-1));
        Assert.That(predictor.ExperimentalElectrophoreticMobility(0), Is.EqualTo(-1));
    }
}
