using System;
using System.Collections.Generic;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

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
            Assert.That(attr.Checksum, Is.EqualTo(checkSum));
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
            Assert.That(attr.Checksum, Is.EqualTo("DEF456"));
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
        [Test]
        public void Constructor_Sets_UniProtSequenceAttributes()
        {
            // Arrange
            var uniProtAttrs = new UniProtSequenceAttributes(
                length: 100,
                mass: 12345,
                checkSum: "CHK123",
                entryModified: new DateTime(2024, 6, 13),
                sequenceVersion: 2,
                isPrecursor: true,
                fragment: UniProtSequenceAttributes.FragmentType.single
            );

            // Act
            var protein = new Protein(
                accession: "P12345",
                sequence: "MPEPTIDESEQ",
                organism: "Homo sapiens",
                isDecoy: false,
                geneNames: new List<Tuple<string, string>>(),
                name: "Test",
                fullName: "Test Protein",
                isContaminant: false,
                sequenceVariations: new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>(),
                spliceSites: new List<SpliceSite>(),
                databaseReferences: new List<DatabaseReference>(),
                databaseFilePath: "db.fasta",
                uniProtSequenceAttributes: uniProtAttrs,
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null
            );

            // Assert
            Assert.That(protein.UniProtSequenceAttributes, Is.EqualTo(uniProtAttrs));
            Assert.That(protein.UniProtSequenceAttributes.Length, Is.EqualTo(100));
            Assert.That(protein.UniProtSequenceAttributes.Mass, Is.EqualTo(12345));
            Assert.That(protein.UniProtSequenceAttributes.Checksum, Is.EqualTo("CHK123"));
            Assert.That(protein.UniProtSequenceAttributes.EntryModified, Is.EqualTo(new DateTime(2024, 6, 13)));
            Assert.That(protein.UniProtSequenceAttributes.SequenceVersion, Is.EqualTo(2));
            Assert.That(protein.UniProtSequenceAttributes.IsPrecursor, Is.True);
            Assert.That(protein.UniProtSequenceAttributes.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.single));
        }

        [Test]
        public void Constructor_Allows_Null_UniProtSequenceAttributes()
        {
            // Act
            var protein = new Protein(
                accession: "P12345",
                sequence: "MPEPTIDESEQ",
                organism: "Homo sapiens",
                isDecoy: false,
                geneNames: new List<Tuple<string, string>>(),
                name: "Test",
                fullName: "Test Protein",
                isContaminant: false,
                sequenceVariations: new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>(),
                spliceSites: new List<SpliceSite>(),
                databaseReferences: new List<DatabaseReference>(),
                databaseFilePath: "db.fasta",
                uniProtSequenceAttributes: null,
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null
            );

            // Assert
            Assert.That(protein.UniProtSequenceAttributes, Is.Not.Null);
            Assert.That(protein.UniProtSequenceAttributes.Length, Is.EqualTo(11));
            Assert.That(protein.UniProtSequenceAttributes.Mass, Is.EqualTo(1275));
        }

        [Test]
        public void UniProtSequenceAttributes_AreAccessible_And_Mutable_IfSet()
        {
            // Arrange
            var uniProtAttrs = new UniProtSequenceAttributes(
                length: 50,
                mass: 5000,
                checkSum: "CHK999",
                entryModified: new DateTime(2020, 1, 1),
                sequenceVersion: 1
            );
            var protein = new Protein(
                accession: "P67890",
                sequence: "MPEPTIDESEQ",
                organism: "Mus musculus",
                isDecoy: false,
                geneNames: new List<Tuple<string, string>>(),
                name: "Test2",
                fullName: "Test Protein 2",
                isContaminant: false,
                sequenceVariations: new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>(),
                spliceSites: new List<SpliceSite>(),
                databaseReferences: new List<DatabaseReference>(),
                databaseFilePath: "db2.fasta",
                uniProtSequenceAttributes: uniProtAttrs,
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null
            );

            // Act
            protein.UniProtSequenceAttributes.UpdateLengthAttribute(60);
            protein.UniProtSequenceAttributes.UpdateMassAttribute(6000);

            // Assert
            Assert.That(protein.UniProtSequenceAttributes.Length, Is.EqualTo(60));
            Assert.That(protein.UniProtSequenceAttributes.Mass, Is.EqualTo(6000));
        }
        [Test]
        public void SequenceAttributes_Default_IsNull()
        {
            // Arrange
            var entry = new ProteinXmlEntry();

            // Assert
            Assert.That(entry.SequenceAttributes, Is.Null);
        }

        [Test]
        public void SequenceAttributes_CanBeSet_AndRetrieved()
        {
            // Arrange
            var entry = new ProteinXmlEntry();
            var attrs = new UniProtSequenceAttributes(
                length: 42,
                mass: 1234,
                checkSum: "CHK",
                entryModified: new DateTime(2024, 6, 13),
                sequenceVersion: 1
            );

            // Act
            entry.SequenceAttributes = attrs;

            // Assert
            Assert.That(entry.SequenceAttributes, Is.EqualTo(attrs));
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(42));
            Assert.That(entry.SequenceAttributes.Mass, Is.EqualTo(1234));
            Assert.That(entry.SequenceAttributes.Checksum, Is.EqualTo("CHK"));
            Assert.That(entry.SequenceAttributes.EntryModified, Is.EqualTo(new DateTime(2024, 6, 13)));
            Assert.That(entry.SequenceAttributes.SequenceVersion, Is.EqualTo(1));
        }

        [Test]
        public void SequenceAttributes_CanBeSetToNull()
        {
            // Arrange
            var entry = new ProteinXmlEntry();
            var attrs = new UniProtSequenceAttributes(10, 100, "A", DateTime.Now, 1);
            entry.SequenceAttributes = attrs;

            // Act
            entry.SequenceAttributes = null;

            // Assert
            Assert.That(entry.SequenceAttributes, Is.Null);
        }

        [Test]
        public void SequenceAttributes_IsMutable_WhenSet()
        {
            // Arrange
            var entry = new ProteinXmlEntry();
            var attrs = new UniProtSequenceAttributes(10, 100, "A", DateTime.Now, 1);
            entry.SequenceAttributes = attrs;

            // Act
            entry.SequenceAttributes.UpdateLengthAttribute(99);
            entry.SequenceAttributes.UpdateMassAttribute(999);

            // Assert
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(99));
            Assert.That(entry.SequenceAttributes.Mass, Is.EqualTo(999));
        }
    }
}
