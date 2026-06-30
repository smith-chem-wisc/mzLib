using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml.Serialization;
using NUnit.Framework;
using pepXML.Generated;

namespace Test.PepXml
{
    /// <summary>
    /// pepXML_v120.cs is a single, very large XSD-generated DTO tree (namespace
    /// pepXML.Generated) of ~80 classes and ~470 properties — exactly the kind of
    /// auto-generated file whose coverage historically tanks because nobody hand-writes
    /// a set/get test for every property. Following the same reflection-driven approach
    /// already used for the mzIdentML DTOs (MzIdentMlReflectionCoverageTests), a single
    /// test walks every concrete generated type and exercises each property's getter and
    /// setter, including the many optional-value "...Specified" flags. A second test runs
    /// a real XmlSerializer round-trip over a populated object graph: constructing the
    /// serializer validates the Xml* attributes across the whole tree, and the
    /// serialize/deserialize cycle exercises the accessors the way real consumers do.
    /// </summary>
    [TestFixture]
    public class PepXmlReflectionCoverageTests
    {
        [Test]
        public void EveryGeneratedType_PropertyGetSet_RoundTrips()
        {
            Assembly assembly = typeof(msms_pipeline_analysis).Assembly;

            Type[] concreteTypes = assembly.GetTypes()
                .Where(t => t.Namespace == "pepXML.Generated"
                            && t.IsClass
                            && !t.IsAbstract
                            && t.GetConstructor(Type.EmptyTypes) != null)
                .ToArray();

            Assert.That(concreteTypes.Length, Is.GreaterThan(0),
                "No instantiable generated types found in pepXML.Generated.");

            foreach (Type type in concreteTypes)
            {
                object instance = Activator.CreateInstance(type);

                foreach (PropertyInfo property in type.GetProperties(BindingFlags.Public | BindingFlags.Instance))
                {
                    if (!property.CanRead || !property.CanWrite || property.GetIndexParameters().Length > 0)
                        continue;

                    object value = CreateValue(property.PropertyType);

                    property.SetValue(instance, value);            // exercises the setter
                    object readBack = property.GetValue(instance); // exercises the getter

                    Assert.That(readBack, Is.EqualTo(value),
                        $"{type.Name}.{property.Name} did not round-trip.");
                }
            }
        }

        /// <summary>
        /// Every generated enum should expose its declared members and survive a
        /// value -> name -> value parse, which is the path the XmlSerializer uses.
        /// </summary>
        [Test]
        public void EveryGeneratedEnum_RoundTripsThroughItsNames()
        {
            Assembly assembly = typeof(msms_pipeline_analysis).Assembly;

            Type[] enums = assembly.GetTypes()
                .Where(t => t.Namespace == "pepXML.Generated" && t.IsEnum)
                .ToArray();

            Assert.That(enums.Length, Is.GreaterThan(0), "No generated enums found.");

            foreach (Type enumType in enums)
            {
                Array values = Enum.GetValues(enumType);
                Assert.That(values.Length, Is.GreaterThan(0), $"{enumType.Name} has no members.");

                foreach (object value in values)
                {
                    string name = Enum.GetName(enumType, value);
                    Assert.That(name, Is.Not.Null);
                    object parsed = Enum.Parse(enumType, name);
                    Assert.That(parsed, Is.EqualTo(value), $"{enumType.Name}.{name} did not round-trip.");
                }
            }
        }

        /// <summary>
        /// A populated msms_pipeline_analysis tree should survive a full XmlSerializer
        /// write/read cycle. Constructing the serializer alone validates the Xml*
        /// attributes on the entire generated tree; the round-trip then confirms the
        /// nested data (run summary -> search summary -> spectrum query -> search hit)
        /// is preserved.
        /// </summary>
        [Test]
        public void MsmsPipelineAnalysis_XmlSerializerRoundTrip_PreservesData()
        {
            var analysis = BuildRepresentativeAnalysis();

            var serializer = new XmlSerializer(typeof(msms_pipeline_analysis));

            string xml;
            using (var writer = new StringWriter())
            {
                serializer.Serialize(writer, analysis);
                xml = writer.ToString();
            }

            Assert.That(xml, Does.Contain("msms_pipeline_analysis"));
            Assert.That(xml, Does.Contain("PEPTIDEK"));

            msms_pipeline_analysis roundTripped;
            using (var reader = new StringReader(xml))
            {
                roundTripped = (msms_pipeline_analysis)serializer.Deserialize(reader);
            }

            Assert.That(roundTripped.name, Is.EqualTo(analysis.name));
            Assert.That(roundTripped.msms_run_summary, Has.Length.EqualTo(1));

            var run = roundTripped.msms_run_summary[0];
            Assert.That(run.base_name, Is.EqualTo("run1"));
            Assert.That(run.sample_enzyme.name, Is.EqualTo("trypsin"));
            Assert.That(run.search_summary, Has.Length.EqualTo(1));
            Assert.That(run.search_summary[0].search_engine, Is.EqualTo(engineType.Comet));

            var query = run.spectrum_query[0];
            Assert.That(query.spectrum, Is.EqualTo("scan=1"));
            Assert.That(query.assumed_charge, Is.EqualTo("2"));

            var hit = query.search_result[0].search_hit[0];
            Assert.That(hit.peptide, Is.EqualTo("PEPTIDEK"));
            Assert.That(hit.hit_rank, Is.EqualTo(1u));
            Assert.That(hit.calc_neutral_pep_mass, Is.EqualTo(927.46f).Within(0.01f));
            Assert.That(hit.modification_info, Is.Not.Null);
            Assert.That(hit.modification_info.mod_aminoacid_mass, Has.Length.EqualTo(1));
            Assert.That(hit.modification_info.mod_aminoacid_mass[0].position, Is.EqualTo("3"));
        }

        private static msms_pipeline_analysis BuildRepresentativeAnalysis()
        {
            return new msms_pipeline_analysis
            {
                name = "test-analysis",
                date = new DateTime(2026, 1, 1),
                summary_xml = "test.pep.xml",
                msms_run_summary = new[]
                {
                    new msms_pipeline_analysisMsms_run_summary
                    {
                        base_name = "run1",
                        raw_data_type = "raw",
                        raw_data = ".mzML",
                        sample_enzyme = new msms_pipeline_analysisMsms_run_summarySample_enzyme
                        {
                            name = "trypsin",
                            specificity = new[]
                            {
                                new msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificity
                                {
                                    cut = "KR",
                                    no_cut = "P",
                                    sense = msms_pipeline_analysisMsms_run_summarySample_enzymeSpecificitySense.C
                                }
                            }
                        },
                        search_summary = new[]
                        {
                            new msms_pipeline_analysisMsms_run_summarySearch_summary
                            {
                                base_name = "run1",
                                search_engine = engineType.Comet,
                                precursor_mass_type = massType.monoisotopic,
                                fragment_mass_type = massType.monoisotopic,
                                search_id = 1,
                                search_database = new msms_pipeline_analysisMsms_run_summarySearch_summarySearch_database
                                {
                                    local_path = "db.fasta",
                                    type = msms_pipeline_analysisMsms_run_summarySearch_summarySearch_databaseType.AA
                                }
                            }
                        },
                        spectrum_query = new[]
                        {
                            new msms_pipeline_analysisMsms_run_summarySpectrum_query
                            {
                                spectrum = "scan=1",
                                start_scan = 1,
                                end_scan = 1,
                                precursor_neutral_mass = 927.46f,
                                assumed_charge = "2",
                                index = 1,
                                search_result = new[]
                                {
                                    new msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_result
                                    {
                                        search_hit = new[]
                                        {
                                            new msms_pipeline_analysisMsms_run_summarySpectrum_querySearch_resultSearch_hit
                                            {
                                                hit_rank = 1,
                                                peptide = "PEPTIDEK",
                                                peptide_prev_aa = "R",
                                                peptide_next_aa = "S",
                                                protein = "PROT1",
                                                num_tot_proteins = 1,
                                                calc_neutral_pep_mass = 927.46f,
                                                massdiff = "0.001",
                                                modification_info = new modInfoDataType
                                                {
                                                    modified_peptide = "PEP[+80]TIDEK",
                                                    mod_aminoacid_mass = new[]
                                                    {
                                                        new modInfoDataTypeMod_aminoacid_mass
                                                        {
                                                            position = "3",
                                                            mass = 181.014
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            };
        }

        /// <summary>
        /// Produces a representative value for a property type so the setter/getter can be
        /// exercised. Abstract/interface reference types that can't be instantiated fall back
        /// to null, which still exercises the accessors.
        /// </summary>
        private static object CreateValue(Type type)
        {
            Type underlying = Nullable.GetUnderlyingType(type);
            if (underlying != null)
                return CreateValue(underlying);

            if (type == typeof(string))
                return "x";
            if (type == typeof(bool))
                return true;
            if (type == typeof(DateTime))
                return new DateTime(2020, 1, 1);
            if (type.IsEnum)
                return Enum.GetValues(type).GetValue(0);
            if (type == typeof(int) || type == typeof(long) || type == typeof(short) || type == typeof(byte) ||
                type == typeof(uint) || type == typeof(ulong) || type == typeof(ushort) || type == typeof(sbyte))
                return Convert.ChangeType(1, type);
            if (type == typeof(double) || type == typeof(float) || type == typeof(decimal))
                return Convert.ChangeType(1, type);
            if (type.IsValueType)
                return Activator.CreateInstance(type); // any other struct -> default
            if (type.IsArray)
                return Array.CreateInstance(type.GetElementType(), 0);
            if (type == typeof(object))
                return "x";

            // Reference type: instantiate when there is an accessible parameterless ctor.
            if (!type.IsAbstract && !type.IsInterface && type.GetConstructor(Type.EmptyTypes) != null)
            {
                try { return Activator.CreateInstance(type); }
                catch { return null; }
            }

            return null; // abstract/interface/no-ctor -> leave null; accessors still run
        }
    }
}
