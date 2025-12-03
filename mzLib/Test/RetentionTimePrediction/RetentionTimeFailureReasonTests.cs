using NUnit.Framework;
using Chromatography.RetentionTimePrediction.Util;
using System.Linq;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for RetentionTimeFailureReason enum
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimeFailureReasonTests
    {
        [Test]
        public void Enum_HasSequenceTooLongValue()
        {
            var reason = RetentionTimeFailureReason.SequenceTooLong;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong));
        }

        [Test]
        public void Enum_HasSequenceTooShortValue()
        {
            var reason = RetentionTimeFailureReason.SequenceTooShort;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooShort));
        }

        [Test]
        public void Enum_HasInvalidAminoAcidValue()
        {
            var reason = RetentionTimeFailureReason.InvalidAminoAcid;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.InvalidAminoAcid));
        }

        [Test]
        public void Enum_HasInvalidMassValue()
        {
            var reason = RetentionTimeFailureReason.InvalidMass;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.InvalidMass));
        }

        [Test]
        public void Enum_HasIncompatibleModificationsValue()
        {
            var reason = RetentionTimeFailureReason.IncompatibleModifications;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));
        }

        [Test]
        public void Enum_HasEmptySequenceValue()
        {
            var reason = RetentionTimeFailureReason.EmptySequence;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
        }

        [Test]
        public void Enum_HasPredictionErrorValue()
        {
            var reason = RetentionTimeFailureReason.PredictionError;
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        [Test]
        public void Enum_CanBeConvertedToString()
        {
            var reason = RetentionTimeFailureReason.SequenceTooLong;
            
            var stringValue = reason.ToString();
            
            Assert.That(stringValue, Is.EqualTo("SequenceTooLong"));
        }

        [Test]
        public void Enum_CanBeParsedFromString()
        {
            var reason = System.Enum.Parse<RetentionTimeFailureReason>("InvalidAminoAcid");
            
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.InvalidAminoAcid));
        }

        [Test]
        public void Enum_ValuesAreDistinct()
        {
            var values = System.Enum.GetValues<RetentionTimeFailureReason>();
            
            Assert.That(values.Length, Is.GreaterThan(0));
            Assert.That(values.Distinct().Count(), Is.EqualTo(values.Length));
        }

        [Test]
        public void Enum_CanBeUsedInSwitch()
        {
            var reason = RetentionTimeFailureReason.IncompatibleModifications;
            string result = reason switch
            {
                RetentionTimeFailureReason.SequenceTooLong => "Too long",
                RetentionTimeFailureReason.SequenceTooShort => "Too short",
                RetentionTimeFailureReason.InvalidAminoAcid => "Invalid AA",
                RetentionTimeFailureReason.InvalidMass => "Invalid mass",
                RetentionTimeFailureReason.IncompatibleModifications => "Incompatible mods",
                RetentionTimeFailureReason.EmptySequence => "Empty",
                RetentionTimeFailureReason.PredictionError => "Error",
                _ => "Unknown"
            };
            
            Assert.That(result, Is.EqualTo("Incompatible mods"));
        }

        [Test]
        public void Enum_NullableCanBeNull()
        {
            RetentionTimeFailureReason? nullableReason = null;
            
            Assert.That(nullableReason, Is.Null);
            Assert.That(nullableReason.HasValue, Is.False);
        }

        [Test]
        public void Enum_NullableCanHaveValue()
        {
            RetentionTimeFailureReason? nullableReason = RetentionTimeFailureReason.SequenceTooLong;
            
            Assert.That(nullableReason, Is.Not.Null);
            Assert.That(nullableReason.HasValue, Is.True);
            Assert.That(nullableReason.Value, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong));
        }

        [Test]
        public void Enum_CanBeCompared()
        {
            var reason1 = RetentionTimeFailureReason.SequenceTooLong;
            var reason2 = RetentionTimeFailureReason.SequenceTooLong;
            var reason3 = RetentionTimeFailureReason.InvalidAminoAcid;
            
            Assert.That(reason1, Is.EqualTo(reason2));
            Assert.That(reason1, Is.Not.EqualTo(reason3));
        }

        [Test]
        public void Enum_AllValuesHaveNames()
        {
            var values = System.Enum.GetValues<RetentionTimeFailureReason>();
            
            foreach (var value in values)
            {
                var name = System.Enum.GetName(value);
                Assert.That(name, Is.Not.Null);
                Assert.That(name, Is.Not.Empty);
            }
        }

        [Test]
        public void Enum_CanBeUsedInOutParameter()
        {
            RetentionTimeFailureReason? reason;
            var success = TryOperation(false, out reason);
            
            Assert.That(success, Is.False);
            Assert.That(reason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        private bool TryOperation(bool shouldSucceed, out RetentionTimeFailureReason? reason)
        {
            if (shouldSucceed)
            {
                reason = null;
                return true;
            }
            else
            {
                reason = RetentionTimeFailureReason.PredictionError;
                return false;
            }
        }
    }
}
