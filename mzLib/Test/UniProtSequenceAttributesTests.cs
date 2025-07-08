using System;
using System.Collections.Generic;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class UniProtSequenceAttributesTests
    {
        [Test]
        public void Constructor_SetsAllMandatoryProperties()
        {
            // Arrange
            int length = 100;
            int mass = 12345;
            string checkSum = "ABC123";
            DateTime entryModified = new DateTime(2024, 6, 13);
            int sequenceVersion = 2;

            // Act
            var attr = new UniProtSequenceAttributes(length, mass, checkSum, entryModified, sequenceVersion);

            // Assert
            Assert.That(attr.Length, Is.EqualTo(length));
            Assert.That(attr.Mass, Is.EqualTo(mass));
            Assert.That(attr.CheckSum, Is.EqualTo(checkSum));
            Assert.That(attr.EntryModified, Is.EqualTo(entryModified));
            Assert.That(attr.SequenceVersion, Is.EqualTo(sequenceVersion));
            Assert.That(attr.IsPrecursor, Is.Null);
            Assert.That(attr.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.unspecified));
        }

        [Test]
        public void Constructor_SetsOptionalProperties()
        {
            // Arrange
            int length = 200;
            int mass = 54321;
            string checkSum = "XYZ789";
            DateTime entryModified = new DateTime(2023, 1, 1);
            int sequenceVersion = 5;
            bool? isPrecursor = true;
            var fragment = UniProtSequenceAttributes.FragmentType.multiple;

            // Act
            var attr = new UniProtSequenceAttributes(length, mass, checkSum, entryModified, sequenceVersion, isPrecursor, fragment);

            // Assert
            Assert.That(attr.IsPrecursor, Is.EqualTo(isPrecursor));
            Assert.That(attr.Fragment, Is.EqualTo(fragment));
        }

        [Test]
        public void Properties_CanBeSetAndGet()
        {
            // Arrange
            var attr = new UniProtSequenceAttributes(10, 20, "DEF456", new DateTime(2022, 12, 31), 4,false,UniProtSequenceAttributes.FragmentType.single);

            // Assert
            Assert.That(attr.Length, Is.EqualTo(10));
            Assert.That(attr.Mass, Is.EqualTo(20));
            Assert.That(attr.CheckSum, Is.EqualTo("DEF456"));
            Assert.That(attr.EntryModified, Is.EqualTo(new DateTime(2022, 12, 31)));
            Assert.That(attr.SequenceVersion, Is.EqualTo(4));
            Assert.That(attr.IsPrecursor, Is.False);
            Assert.That(attr.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.single));
        }

        [Test]
        public void FragmentType_Enum_Values_AreDistinct()
        {
            Assert.That(UniProtSequenceAttributes.FragmentType.single, Is.Not.EqualTo(UniProtSequenceAttributes.FragmentType.multiple));
            Assert.That(UniProtSequenceAttributes.FragmentType.unspecified, Is.Not.EqualTo(UniProtSequenceAttributes.FragmentType.single));
            Assert.That(UniProtSequenceAttributes.FragmentType.unspecified, Is.Not.EqualTo(UniProtSequenceAttributes.FragmentType.multiple));
        }
        [Test]
        public void UpdateLengthAttribute_WithInt_UpdatesLength()
        {
            // Arrange
            var attr = new UniProtSequenceAttributes(10, 1000, "CHK", DateTime.Now, 1);

            // Act
            attr.UpdateLengthAttribute(25);

            // Assert
            Assert.That(attr.Length, Is.EqualTo(25));
        }

        [Test]
        public void UpdateLengthAttribute_WithString_UpdatesLengthToStringLength()
        {
            // Arrange
            var attr = new UniProtSequenceAttributes(10, 1000, "CHK", DateTime.Now, 1);

            // Act
            attr.UpdateLengthAttribute("MPEPTIDESEQ");

            // Assert
            Assert.That(attr.Length, Is.EqualTo("MPEPTIDESEQ".Length));
        }
        [Test]
        public void UpdateMassAttribute_WithInt_UpdatesMass()
        {
            // Arrange
            var attr = new UniProtSequenceAttributes(10, 1000, "CHK", DateTime.Now, 1);

            // Act
            attr.UpdateMassAttribute(2500);

            // Assert
            Assert.That(attr.Mass, Is.EqualTo(2500));
        }

        [Test]
        public void UpdateMassAttribute_WithString_UpdatesMassToMonoisotopicMass()
        {
            // Arrange
            var attr = new UniProtSequenceAttributes(10, 1000, "CHK", DateTime.Now, 1);
            string sequence = "PEPTIDE";
            // The expected mass is calculated using PeptideWithSetModifications
            var peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());
            int expectedMass = (int)Math.Round(peptide.MonoisotopicMass);

            // Act
            attr.UpdateMassAttribute(sequence);

            // Assert
            Assert.That(attr.Mass, Is.EqualTo(expectedMass));
        }
    }
}
