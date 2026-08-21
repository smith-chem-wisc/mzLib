using System;
using System.IO;
using System.Xml;
using System.Xml.Serialization;
using NUnit.Framework;

namespace Test.MzIdentML
{
    /// <summary>
    /// The mzIdentML schema version and its target namespace are not the same string, which is easy
    /// to get wrong: 1.1.0 declares targetNamespace .../mzIdentML/1.1, while 1.1.1 declares
    /// .../1.1.1 and 1.2.0 declares .../1.2 (HUPO-PSI/mzIdentML, schema/*.xsd).
    ///
    /// Getting it wrong on the root type is silent in both directions -- a wrong namespace still
    /// serializes and still round-trips through our own deserializer, it just makes the file
    /// something no other mzIdentML consumer recognises, and makes standards-compliant files
    /// unreadable here.
    /// </summary>
    [TestFixture]
    public sealed class MzIdentMlNamespaceTests
    {
        private const string Namespace110 = "http://psidev.info/psi/pi/mzIdentML/1.1";
        private const string Namespace111 = "http://psidev.info/psi/pi/mzIdentML/1.1.1";
        private const string Namespace120 = "http://psidev.info/psi/pi/mzIdentML/1.2";

        [Test]
        [TestCase(typeof(mzIdentML110.Generated.MzIdentMLType110), Namespace110)]
        [TestCase(typeof(mzIdentML111.Generated.MzIdentMLType111), Namespace111)]
        [TestCase(typeof(mzIdentML120.Generated.MzIdentMLType120), Namespace120)]
        public void RootTypeDeclaresTheSchemaTargetNamespace(Type rootType, string expectedNamespace)
        {
            var root = (XmlRootAttribute)Attribute.GetCustomAttribute(rootType, typeof(XmlRootAttribute));
            var type = (XmlTypeAttribute)Attribute.GetCustomAttribute(rootType, typeof(XmlTypeAttribute));

            Assert.That(root, Is.Not.Null, $"{rootType.Name} has no XmlRootAttribute");
            Assert.That(root.ElementName, Is.EqualTo("MzIdentML"));
            Assert.That(root.Namespace, Is.EqualTo(expectedNamespace));
            Assert.That(type?.Namespace, Is.EqualTo(expectedNamespace));
        }

        [Test]
        [TestCase(typeof(mzIdentML110.Generated.MzIdentMLType110), Namespace110)]
        [TestCase(typeof(mzIdentML111.Generated.MzIdentMLType111), Namespace111)]
        [TestCase(typeof(mzIdentML120.Generated.MzIdentMLType120), Namespace120)]
        public void SerializedRootElementIsInTheSchemaTargetNamespace(Type rootType, string expectedNamespace)
        {
            object instance = Activator.CreateInstance(rootType);

            using var buffer = new MemoryStream();
            new XmlSerializer(rootType).Serialize(buffer, instance);
            buffer.Position = 0;

            using var reader = XmlReader.Create(buffer);
            while (reader.Read())
            {
                if (reader.NodeType != XmlNodeType.Element) continue;

                Assert.That(reader.LocalName, Is.EqualTo("MzIdentML"));
                Assert.That(reader.NamespaceURI, Is.EqualTo(expectedNamespace));
                return;
            }

            Assert.Fail("no root element was serialized");
        }
    }
}
