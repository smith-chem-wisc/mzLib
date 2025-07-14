using NUnit.Framework;
using System;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class ProteinDbLoaderTests
    {
        [Test]
        public void SanitizeAminoAcidSequence_ValidSequence_ReturnsUnchanged()
        {
            string input = "ACDEFGHIKLMNPQRSTVWY";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(input));
        }

        [Test]
        public void SanitizeAminoAcidSequence_SequenceWithUnicode_ReplacesWithX()
        {
            string input = "ACD\u03A9EFG"; // \u03A9 is Greek Omega
            string expected = "ACDXEFG";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void SanitizeAminoAcidSequence_SequenceWithInvalidCharacters_ReplacesWithX()
        {
            string input = "ACD*EFG#HIK";
            string expected = "ACDXEFGXHIK";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void SanitizeAminoAcidSequence_SequenceWithWhitespace_ReplacesWithX()
        {
            string input = "ACD EFG\tHIK";
            string expected = "ACDXEFGXHIK";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void SanitizeAminoAcidSequence_EmptyString_ReturnsEmpty()
        {
            string input = "";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(""));
        }

        [Test]
        public void SanitizeAminoAcidSequence_AllInvalidCharacters_ReplacesAllWithX()
        {
            string input = "!@#$%^&*()";
            string expected = new string('X', input.Length);
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void SanitizeAminoAcidSequence_NullInput_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() =>
            {
                ProteinDbLoader.SanitizeAminoAcidSequence(null, 'X');
            });
        }

        [Test]
        public void SanitizeAminoAcidSequence_ValidAndInvalidMix_ReplacesOnlyInvalid()
        {
            string input = "ACD*EFGZ";
            // Assuming Z is not a valid residue
            string expected = "ACDXEFGX";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void SanitizeAminoAcidSequence_ValidSequenceWithLowercase_ReplacesLowercaseWithX()
        {
            string input = "ACDefgHIK";
            string expected = "ACDXXXHIK";
            string result = ProteinDbLoader.SanitizeAminoAcidSequence(input, 'X');
            Assert.That(result, Is.EqualTo(expected));
        }
    }
}
