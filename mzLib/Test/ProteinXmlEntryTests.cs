using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class ProteinXmlEntryTests
    {
        private static XmlReader CreateSequenceReader(string attributes, string sequence = "ACDE")
        {
            string xml = $"<sequence {attributes}>{sequence}</sequence>";
            return XmlReader.Create(new StringReader(xml));
        }

        [Test]
        public void ParseSequenceAttributes_ParsesChecksum()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("checksum=\"ABC123\"");
            reader.Read(); // Move to <sequence>
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That("ABC123", Is.EqualTo(entry.SequenceAttributes.Checksum));
        }

        [Test]
        public void ParseSequenceAttributes_ParsesModifiedDate()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("modified=\"2022-01-15\"");
            reader.Read();
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That(new DateTime(2022, 1, 15), Is.EqualTo(entry.SequenceAttributes.EntryModified));
        }

        [Test]
        public void ParseSequenceAttributes_ParsesVersion()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("version=\"7\"");
            reader.Read();
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That(7, Is.EqualTo(entry.SequenceAttributes.SequenceVersion));
        }

        [Test]
        public void ParseSequenceAttributes_ParsesPrecursor()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("precursor=\"true\"");
            reader.Read();
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That(entry.SequenceAttributes.IsPrecursor.Value);
        }

        [Test]
        public void ParseSequenceAttributes_ParsesFragment()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("fragment=\"single\"");
            reader.Read();
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That(UniProtSequenceAttributes.FragmentType.single, Is.EqualTo(entry.SequenceAttributes.Fragment));
        }

        [Test]
        public void ParseSequenceAttributes_ComputesLengthAndMass()
        {
            var entry = new ProteinXmlEntry();
            using var reader = CreateSequenceReader("");
            reader.Read();
            entry.GetType().GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(entry, new object[] { reader });
            Assert.That(4, Is.EqualTo(entry.SequenceAttributes.Length));
            Assert.That(entry.SequenceAttributes.Mass, Is.GreaterThan(0));
        }
        private ProteinXmlEntry CreateEntryAndParse(string xml)
        {
            var entry = new ProteinXmlEntry();
            using var reader = XmlReader.Create(new StringReader(xml));
            reader.Read(); // Move to <sequence>
            var method = typeof(ProteinXmlEntry).GetMethod("ParseSequenceAttributes", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            method.Invoke(entry, new object[] { reader });
            return entry;
        }

        [Test]
        public void Parses_All_Attributes_Correctly()
        {
            string xml = "<sequence length=\"10\" mass=\"1234\" checksum=\"CHK\" modified=\"2024-06-13\" version=\"2\" precursor=\"true\" fragment=\"single\">PEPTIDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(7)); //note that the read length 10 is incorrect for the sequence PEPTIDE. It should be corrected to 7
            Assert.That(entry.SequenceAttributes.Mass, Is.EqualTo(799)); //note that the read mass 1234 is incorrect for the sequence PEPTIDE. It should be corrected to 799
            Assert.That(entry.SequenceAttributes.Checksum, Is.EqualTo("CHK"));
            Assert.That(entry.SequenceAttributes.EntryModified, Is.EqualTo(new DateTime(2024, 6, 13)));
            Assert.That(entry.SequenceAttributes.SequenceVersion, Is.EqualTo(2));
            Assert.That(entry.SequenceAttributes.IsPrecursor, Is.True);
            Assert.That(entry.SequenceAttributes.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.single));
        }

        [Test]
        public void Handles_Missing_Optional_Attributes()
        {
            string xml = "<sequence>ACDEFGHIKL</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(10));
            Assert.That(entry.Sequence, Is.EqualTo("ACDEFGHIKL"));
            Assert.That(entry.SequenceAttributes.IsPrecursor, Is.False);
            Assert.That(entry.SequenceAttributes.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.unspecified));
        }

        [Test]
        public void Handles_Malformed_Length_And_Mass()
        {
            string xml = "<sequence length=\"abc\" mass=\"xyz\">ACDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(4)); // fallback to sequence length
            Assert.That(entry.SequenceAttributes.Mass, Is.GreaterThan(0)); // fallback to computed mass
        }

        [Test]
        public void Handles_Missing_Sequence_Text()
        {
            string xml = "<sequence length=\"5\" mass=\"1000\" />";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.Length, Is.EqualTo(0)); //This value is 0 'zero' because no sequence text is provided and it is computed rather than read
            Assert.That(entry.SequenceAttributes.Mass, Is.EqualTo(0)); //This value is 0 'zero' because no sequence text is provided and it is computed rather than read
            Assert.That(entry.Sequence, Is.Empty);
        }

        [Test]
        public void Handles_Invalid_Modified_Date()
        {
            string xml = "<sequence modified=\"notadate\">ACDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);

            // Compare only the date part (year, month, day)
            var expected = DateTime.Now.Date;
            var actual = entry.SequenceAttributes.EntryModified.Date;
            Assert.That(actual, Is.EqualTo(expected));
        }

        [Test]
        public void Handles_Invalid_Version()
        {
            string xml = "<sequence version=\"notanint\">ACDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.SequenceVersion, Is.EqualTo(-1));
        }

        [Test]
        public void Handles_Invalid_Precursor()
        {
            string xml = "<sequence precursor=\"notabool\">ACDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.IsPrecursor, Is.False);
        }

        [Test]
        public void Handles_Invalid_Fragment()
        {
            string xml = "<sequence fragment=\"notavalidfragment\">ACDE</sequence>";
            var entry = CreateEntryAndParse(xml);
            Assert.That(entry.SequenceAttributes, Is.Not.Null);
            Assert.That(entry.SequenceAttributes.Fragment, Is.EqualTo(UniProtSequenceAttributes.FragmentType.unspecified));
        }
    }
}
