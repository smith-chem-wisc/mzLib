using NUnit.Framework;
using Chromatography.RetentionTimePrediction.Util;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for IncompatibleModHandlingMode enum
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class IncompatibleModHandlingModeTests
    {
        [Test]
        public void Enum_HasRemoveIncompatibleModsValue()
        {
            var mode = IncompatibleModHandlingMode.RemoveIncompatibleMods;
            
            Assert.That(mode, Is.EqualTo(IncompatibleModHandlingMode.RemoveIncompatibleMods));
        }

        [Test]
        public void Enum_HasUsePrimarySequenceValue()
        {
            var mode = IncompatibleModHandlingMode.UsePrimarySequence;
            
            Assert.That(mode, Is.EqualTo(IncompatibleModHandlingMode.UsePrimarySequence));
        }

        [Test]
        public void Enum_HasThrowExceptionValue()
        {
            var mode = IncompatibleModHandlingMode.ThrowException;
            
            Assert.That(mode, Is.EqualTo(IncompatibleModHandlingMode.ThrowException));
        }

        [Test]
        public void Enum_HasReturnNullValue()
        {
            var mode = IncompatibleModHandlingMode.ReturnNull;
            
            Assert.That(mode, Is.EqualTo(IncompatibleModHandlingMode.ReturnNull));
        }

        [Test]
        public void Enum_CanBeConvertedToString()
        {
            var mode = IncompatibleModHandlingMode.RemoveIncompatibleMods;
            
            var stringValue = mode.ToString();
            
            Assert.That(stringValue, Is.EqualTo("RemoveIncompatibleMods"));
        }

        [Test]
        public void Enum_CanBeParsedFromString()
        {
            var mode = System.Enum.Parse<IncompatibleModHandlingMode>("UsePrimarySequence");
            
            Assert.That(mode, Is.EqualTo(IncompatibleModHandlingMode.UsePrimarySequence));
        }

        [Test]
        public void Enum_CanBeUsedInSwitch()
        {
            var mode = IncompatibleModHandlingMode.ThrowException;
            string result = mode switch
            {
                IncompatibleModHandlingMode.RemoveIncompatibleMods => "Remove",
                IncompatibleModHandlingMode.UsePrimarySequence => "Primary",
                IncompatibleModHandlingMode.ThrowException => "Throw",
                IncompatibleModHandlingMode.ReturnNull => "Null",
                _ => "Unknown"
            };
            
            Assert.That(result, Is.EqualTo("Throw"));
        }

        [Test]
        public void Enum_CanBeCompared()
        {
            var mode1 = IncompatibleModHandlingMode.RemoveIncompatibleMods;
            var mode2 = IncompatibleModHandlingMode.RemoveIncompatibleMods;
            var mode3 = IncompatibleModHandlingMode.ThrowException;
            
            Assert.That(mode1, Is.EqualTo(mode2));
            Assert.That(mode1, Is.Not.EqualTo(mode3));
        }

        [Test]
        public void Enum_AllValuesHaveNames()
        {
            var values = System.Enum.GetValues<IncompatibleModHandlingMode>();
            
            foreach (var value in values)
            {
                var name = System.Enum.GetName(value);
                Assert.That(name, Is.Not.Null);
                Assert.That(name, Is.Not.Empty);
            }
        }
    }
}
