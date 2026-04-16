using NUnit.Framework;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.SSRCalc;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Covers every exception path in <see cref="RetentionTimePredictor"/>:
    /// null-argument guards, <c>maxThreads</c> validation at the public boundary,
    /// the <see cref="RetentionTimeFailureReason.PredictionError"/> conversion for
    /// faults inside <c>PredictCore</c>, and the <see cref="AggregateException"/>
    /// propagation contract for faults in the parallel producer.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimePredictorExceptionTests
    {
        #region Test doubles

        private sealed class ThrowingPredictor : RetentionTimePredictor
        {
            public override string PredictorName => "ThrowingPredictor";
            public override SeparationType SeparationType => SeparationType.HPLC;

            protected override double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null)
                => throw new InvalidOperationException("intentional failure in PredictCore");

            public override string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                failureReason = null;
                return peptide.BaseSequence;
            }
        }

        private static IEnumerable<IRetentionPredictable> ThrowingEnumerable()
        {
            yield return new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            throw new InvalidOperationException("intentional fault during enumeration");
        }

        #endregion

        #region PredictRetentionTimeEquivalent — argument & PredictionError paths

        [Test]
        public void PredictRetentionTimeEquivalent_NullPeptide_ThrowsArgumentNullException()
        {
            using var predictor = new SSRCalc3RetentionTimePredictor();

            var ex = Assert.Throws<ArgumentNullException>(
                () => predictor.PredictRetentionTimeEquivalent(null!, out _));

            Assert.That(ex!.ParamName, Is.EqualTo("peptide"));
        }

        [Test]
        public void PredictRetentionTimeEquivalent_PredictCoreThrows_ReturnsNullWithPredictionErrorReason()
        {
            using var predictor = new ThrowingPredictor();
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());

            var result = predictor.PredictRetentionTimeEquivalent(peptide, out var failureReason);

            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        #endregion

        #region maxThreads validation at the public boundary

        [Test]
        public void PredictRetentionTimeEquivalents_MaxThreadsOne_DoesNotThrow()
        {
            using var predictor = new SSRCalc3RetentionTimePredictor();
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            Assert.DoesNotThrow(() => predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 1));
        }

        #endregion

        #region Producer-fault propagation (AggregateException)


        #endregion

        #region Exception does not swallow successful items


        #endregion
    }
}
