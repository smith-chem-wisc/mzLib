using NUnit.Framework;
using Chromatography;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for SeparationType enum
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SeparationTypeTests
    {
        [Test]
        public void Enum_HasHPLCValue()
        {
            var type = SeparationType.HPLC;
            
            Assert.That(type, Is.EqualTo(SeparationType.HPLC));
        }

        [Test]
        public void Enum_HasCZEValue()
        {
            var type = SeparationType.CZE;
            
            Assert.That(type, Is.EqualTo(SeparationType.CZE));
        }

        [Test]
        public void Enum_CanBeConvertedToString()
        {
            var type = SeparationType.HPLC;
            
            var stringValue = type.ToString();
            
            Assert.That(stringValue, Is.EqualTo("HPLC"));
        }

        [Test]
        public void Enum_CanBeParsedFromString()
        {
            var type = System.Enum.Parse<SeparationType>("CZE");
            
            Assert.That(type, Is.EqualTo(SeparationType.CZE));
        }

        [Test]
        public void Enum_CanBeUsedInSwitch()
        {
            var type = SeparationType.HPLC;
            string result = type switch
            {
                SeparationType.HPLC => "High Performance Liquid Chromatography",
                SeparationType.CZE => "Capillary Zone Electrophoresis",
                _ => "Unknown"
            };
            
            Assert.That(result, Is.EqualTo("High Performance Liquid Chromatography"));
        }

        [Test]
        public void Enum_CanBeCompared()
        {
            var type1 = SeparationType.HPLC;
            var type2 = SeparationType.HPLC;
            var type3 = SeparationType.CZE;
            
            Assert.That(type1, Is.EqualTo(type2));
            Assert.That(type1, Is.Not.EqualTo(type3));
        }

        [Test]
        public void Enum_AllValuesHaveNames()
        {
            var values = System.Enum.GetValues<SeparationType>();
            
            foreach (var value in values)
            {
                var name = System.Enum.GetName(value);
                Assert.That(name, Is.Not.Null);
                Assert.That(name, Is.Not.Empty);
            }
        }
    }
}
