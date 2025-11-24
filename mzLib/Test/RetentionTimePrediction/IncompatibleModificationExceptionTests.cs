using NUnit.Framework;
using Chromatography.RetentionTimePrediction.Util;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for IncompatibleModificationException
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class IncompatibleModificationExceptionTests
    {
        [Test]
        public void Constructor_WithSequenceAndPredictorName_SetsProperties()
        {
            var exception = new IncompatibleModificationException(
                "PEPTIDE[HexNAc]",
                "PEPTIDE[+203.07]",
                "TestPredictor");
            
            Assert.That(exception.OriginalSequence, Is.EqualTo("PEPTIDE[HexNAc]"));
            Assert.That(exception.WorkingSequence, Is.EqualTo("PEPTIDE[+203.07]"));
            Assert.That(exception.PredictorName, Is.EqualTo("TestPredictor"));
        }

        [Test]
        public void Message_ContainsRelevantInformation()
        {
            var exception = new IncompatibleModificationException(
                "PEPTIDE[HexNAc]",
                "PEPTIDE[+203.07]",
                "Chronologer");
            
            var message = exception.Message;
            
            Assert.That(message, Does.Contain("Chronologer"));
            Assert.That(message, Does.Contain("PEPTIDE"));
        }

        [Test]
        public void Constructor_EmptyStrings_DoesNotThrow()
        {
            Assert.DoesNotThrow(() => new IncompatibleModificationException("", "", ""));
        }

        [Test]
        public void Constructor_NullStrings_DoesNotThrow()
        {
            Assert.DoesNotThrow(() => new IncompatibleModificationException(null, null, null));
        }

        [Test]
        public void Exception_CanBeCaught_AsBaseException()
        {
            try
            {
                throw new IncompatibleModificationException("SEQ", "SEQ", "Pred");
            }
            catch (System.Exception ex)
            {
                Assert.That(ex, Is.InstanceOf<IncompatibleModificationException>());
            }
        }

        [Test]
        public void Exception_PreservesStackTrace()
        {
            try
            {
                ThrowException();
                Assert.Fail("Should have thrown exception");
            }
            catch (IncompatibleModificationException ex)
            {
                Assert.That(ex.StackTrace, Is.Not.Null);
                Assert.That(ex.StackTrace, Does.Contain(nameof(ThrowException)));
            }
        }

        private void ThrowException()
        {
            throw new IncompatibleModificationException("SEQ", "SEQ", "Pred");
        }
    }
}
