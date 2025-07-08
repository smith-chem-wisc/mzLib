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
    }
}
