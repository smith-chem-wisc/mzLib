using NUnit.Framework;
using System;
using System.IO;
using System.Text;
using System.Xml;
using System.Xml.Serialization;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestMzIdentML
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        //[Test]
        //[TestCase("MPC_1_1_example_Multiple_search_engines.mzid", "1.1.0")]
        //[TestCase("MPC_1_1_example_Multiple_search_engines.mzid", "1.1.1")]
        //[TestCase("combined_1.2.mzid", "1.2.0")]
        //public static void ReadMzIDTest(string filename, string version)
        //{
        //    string mzidfile = Path.Combine(Environment.CurrentDirectory, "DatabaseTests", "mzIdentML", filename);
        //    Type providerType = version.StartsWith("1.1") ? typeof(mzIdentML110.Generated.ProviderType) : typeof(mzIdentML120.Generated.ProviderType);
        //    var mySerializer = new XmlSerializer(providerType);
        //    using var myFileStream = new FileStream(mzidfile, FileMode.Open);
        //    if (version == "1.1.0")
        //    {
        //        var myObject = (mzIdentML110.Generated.ProviderType)mySerializer.Deserialize(myFileStream);
        //    }
        //    else if (version == "1.1.1")
        //    {
        //        var myObject = (mzIdentML111.Generated.ProviderType)mySerializer.Deserialize(myFileStream);
        //    }
        //    else if (version == "1.2.0")
        //    {
        //        var myObject = (mzIdentML120.Generated.ProviderType)mySerializer.Deserialize(myFileStream);
        //    }
        //    else
        //    {
        //        throw new ArgumentException("mzID version not supported");
        //    }
        //}

        [Test]
        public static void WriteMzID110Test()
        {
            string version = "1.1.0";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML110.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        [Test]
        public static void WriteMzID111Test()
        {
            string version = "1.1.1";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML111.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML111.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML111.Generated.RoleType()
                    {
                        cvParam = new mzIdentML111.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        [Test]
        public static void WriteMzID120Test()
        {
            string version = "1.2.0";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML120.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML120.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML120.Generated.RoleType()
                    {
                        cvParam = new mzIdentML120.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }
    }
}