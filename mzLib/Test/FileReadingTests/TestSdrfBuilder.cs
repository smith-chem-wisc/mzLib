using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for building an SDRF document from mzLib-native inputs.
    ///
    /// The behaviour these pin hardest is what the builder REFUSES to do: invent a sample fact,
    /// copy an accession from an input, or write a plausible neighbour for a term it cannot
    /// resolve. Each of those produces a document that passes validation and says nothing true.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfBuilder
    {
        private static CvParam Psi(string accession)
        {
            Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession(accession, out var t), Is.True);
            return t;
        }

        private static Modification Mod(string id, string motif, string unimod)
        {
            ModificationMotif.TryGetMotif(motif, out var target);
            return new Modification(
                _originalId: id, _modificationType: "Common", _target: target,
                _locationRestriction: "Anywhere.",
                _databaseReference: new Dictionary<string, IList<string>> { { "Unimod", new List<string> { unimod } } });
        }

        private static SdrfSample Sample(string name = "Sample 1") => new()
        {
            SourceName = name,
            Organism = new CvParam("NCBITaxon", "NCBITaxon:9606", "Homo sapiens", ""),
            Characteristics = new Dictionary<string, CvParam>
            {
                ["characteristics[organism part]"] = new CvParam("UBERON", "UBERON:0002107", "liver", "")
            },
            BiologicalReplicate = 1,
            Label = new CvParam("", "", "label free sample", "")
        };

        private static SdrfAssay Assay(string file = "a.raw") => new()
        {
            DataFileName = file,
            AssayName = "run " + file,
            Instrument = Psi("MS:1001911"),
            PrecursorMassTolerance = new PpmTolerance(5),
            ProductMassTolerance = new PpmTolerance(20),
            CleavageAgent = ProteaseDictionary.Dictionary["trypsin"],
            FixedModifications = new[] { Mod("Carbamidomethyl on C", "C", "4") },
            VariableModifications = new[] { Mod("Oxidation on M", "M", "35") },
            DissociationType = DissociationType.HCD,
            AcquisitionMethod = new CvParam("PRIDE", "PRIDE:0000627", "Data-dependent acquisition", ""),
            TechnicalReplicate = 1,
            Fraction = 1
        };

        [Test]
        public void BuildsAValidDocument()
        {
            var document = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) });
            var result = SdrfValidator.Validate(document);

            Assert.That(result.IsValid, Is.True,
                "errors: " + string.Join(" | ", result.Errors.Select(e => e.ToString())));
            Assert.That(document.Results.Count, Is.EqualTo(1));
        }

        [Test]
        public void WritesEveryColumnAsAControlledVocabularyTerm()
        {
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) }).Results[0];

            Assert.That(row["characteristics[organism]"], Is.EqualTo("NT=Homo sapiens;AC=NCBITaxon:9606"));
            Assert.That(row["comment[instrument]"], Is.EqualTo("NT=Q Exactive;AC=MS:1001911"));
            // MS:1001313 is "Trypsin/P", NOT plain Trypsin (MS:1001251). The builder is faithfully
            // emitting what mzLib's embedded proteases.tsv says, and that row is wrong -- the same
            // bug that makes every mzIdentML MetaMorpheus exports declare Trypsin/P for ordinary
            // trypsin. Pinned deliberately so this test FAILS when the proteases.tsv fix lands,
            // rather than quietly agreeing with whatever the data says.
            Assert.That(row["comment[cleavage agent details]"], Does.Contain("AC=MS:1001313"),
                "expected to change to MS:1001251 once the proteases.tsv trypsin row is corrected");
            Assert.That(row["comment[proteomics data acquisition method]"],
                Is.EqualTo("NT=Data-dependent acquisition;AC=PRIDE:0000627"));
        }

        [Test]
        public void ModificationsCarryAccessionTargetAndFixedOrVariable()
        {
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) }).Results[0];
            var mods = row.All("comment[modification parameters]");

            Assert.That(mods.Count, Is.EqualTo(2));

            // EXACT, not Does.Contain. The contains-form passed for as long as NT= carried
            // "Carbamidomethyl on C" — the motif duplicated into the name field, where the
            // specification puts only the bare name and leaves the residue to TA=. Every other
            // column in this fixture is pinned exactly; these two were the only ones that were not,
            // which is precisely why a green suite could not see it.
            Assert.That(mods[0], Is.EqualTo("NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"));
            Assert.That(mods[1], Is.EqualTo("NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable"));
        }

        /// <summary>
        /// mzLib's Modification carries the motif in IdWithMotif ("Oxidation on M") and the bare
        /// name in OriginalId (Modification.cs:57-84). SDRF wants the bare name in NT= and the
        /// residue in TA=, so a cell named from IdWithMotif states the residue twice and states it
        /// once in a field that is not for it.
        /// </summary>
        [Test]
        public void TheModificationNameIsBareAndTheResidueAppearsOnlyInTa()
        {
            var mod = Mod("Oxidation on M", "M", "35");
            Assert.That(mod.IdWithMotif, Is.EqualTo("Oxidation on M"),
                "the input really does carry the motif; that is what makes this worth pinning");
            Assert.That(mod.OriginalId, Is.EqualTo("Oxidation"));

            var assay = Assay() with
            {
                FixedModifications = Array.Empty<Modification>(),
                VariableModifications = new[] { mod }
            };
            string cell = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) })
                .Results[0]["comment[modification parameters]"];

            Assert.That(cell, Is.EqualTo("NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable"));
            Assert.That(cell, Does.Not.Contain(" on "), "the motif belongs to TA=, not to NT=");
        }

        [Test]
        public void DissociationResolvesFromThePinnedVocabulary()
        {
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) }).Results[0];
            Assert.That(row["comment[dissociation method]"], Does.Contain("AC=MS:1000422"));
        }

        [Test]
        public void TolerancesAreWrittenAsValueAndUnit()
        {
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) }).Results[0];
            Assert.That(row["comment[precursor mass tolerance]"], Is.EqualTo("5 ppm"));
            Assert.That(row["comment[fragment mass tolerance]"], Is.EqualTo("20 ppm"));

            var absolute = Assay() with { ProductMassTolerance = new AbsoluteTolerance(0.02) };
            var daltons = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), absolute) }).Results[0];
            Assert.That(daltons["comment[fragment mass tolerance]"], Is.EqualTo("0.02 Da"));
        }

        /// <summary>
        /// A tolerance is a number in a file format, not a number shown to a person. Under the
        /// machine's own culture this wrote "0,02 Da" on a de-DE box, which no SDRF consumer parses
        /// — and the rest of the suite could not see it, because the machine it runs on is en-US.
        /// </summary>
        [Test]
        [SetCulture("de-DE")]
        public void TolerancesAreWrittenInInvariantCultureWhateverTheMachineIsSetTo()
        {
            Assert.That(0.02d.ToString(), Is.EqualTo("0,02"),
                "guard on the fixture itself: this test proves nothing unless the culture took");

            var assay = Assay() with
            {
                PrecursorMassTolerance = new PpmTolerance(5.5),
                ProductMassTolerance = new AbsoluteTolerance(0.02)
            };
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }).Results[0];

            Assert.That(row["comment[precursor mass tolerance]"], Is.EqualTo("5.5 ppm"));
            Assert.That(row["comment[fragment mass tolerance]"], Is.EqualTo("0.02 Da"));
        }

        /// <summary>
        /// The specification asks for vMAJOR.MINOR.PATCH. Stripping only a lowercase prefix turned a
        /// caller's "V1.1.0" into "vV1.1.0", and nothing downstream validates this column, so the
        /// malformed value went out silently.
        /// </summary>
        [Test]
        public void TheSdrfVersionPrefixIsWrittenOnceWhateverCasingTheCallerSupplied()
        {
            string Version(string supplied) => SdrfBuilder
                .Build(new[] { new SdrfRowInput(Sample(), Assay()) },
                       new SdrfBuilderOptions { SdrfVersion = supplied })
                .Results[0]["comment[sdrf version]"];

            Assert.That(Version("1.1.0"), Is.EqualTo("v1.1.0"));
            Assert.That(Version("v1.1.0"), Is.EqualTo("v1.1.0"));
            Assert.That(Version("V1.1.0"), Is.EqualTo("v1.1.0"));
        }

        // ---------------- what it refuses to do ----------------

        /// <summary>
        /// The decision this whole feature turns on. A builder that pads a missing sample fact with
        /// a reserved word produces a document that passes every validator, generates no drift
        /// findings, and cannot answer a question.
        /// </summary>
        [Test]
        public void AMissingSampleFactFailsRatherThanBecomingNotAvailable()
        {
            var sample = Sample() with { Organism = null };

            var thrown = Assert.Throws<MzLibException>(
                () => SdrfBuilder.Build(new[] { new SdrfRowInput(sample, Assay()) }));

            Assert.That(thrown.Message, Does.Contain("characteristics[organism]"));
            Assert.That(thrown.Message, Does.Contain("cannot be mined"));
        }

        [Test]
        public void OptingOutExplicitlyIsAllowedAndCoverageThenFlagsIt()
        {
            var sample = Sample() with { Organism = null };
            var options = new SdrfBuilderOptions { RequireSampleMetadata = false };

            var document = SdrfBuilder.Build(new[] { new SdrfRowInput(sample, Assay()) }, options);

            Assert.That(document.Results[0]["characteristics[organism]"], Is.EqualTo("not available"));
            Assert.That(SdrfValidator.Validate(document).IsValid, Is.True, "still passes validation");

            var collection = new SdrfCollection(new[] { document }, new[] { "opted-out" });
            Assert.That(SdrfCoverage.Uninformative(collection).Select(c => c.Column),
                Does.Contain("characteristics[organism]"),
                "which is exactly why the coverage metric exists");
        }

        [Test]
        public void AnUnmappableDissociationTypeIsOmittedNotApproximated()
        {
            // aEPD has no PSI-MS term. Writing a near neighbour would silently join this run to the
            // wrong population.
            var assay = Assay() with { DissociationType = DissociationType.aEPD };

            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }),
                Throws.TypeOf<MzLibException>());

            var lenient = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) },
                new SdrfBuilderOptions { RequireSampleMetadata = false });
            Assert.That(lenient.Results[0]["comment[dissociation method]"], Is.EqualTo("not available"));
        }

        [Test]
        public void ZeroBasedReplicatesAreRejected()
        {
            // mzLib's SpectraFileInfo stores these 0-based; SDRF is 1-based. Forwarding directly
            // would write a 0 that every consumer reads as an error.
            var assay = Assay() with { Fraction = 0 };
            var thrown = Assert.Throws<MzLibException>(
                () => SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }));
            Assert.That(thrown.Message, Does.Contain("1-based"));
        }

        [Test]
        public void AThermoNameOnlyInstrumentIsResolvedAgainstPsiMs()
        {
            // A RAW file yields a name with an empty accession. Resolving it here is correct;
            // putting an ontology in the raw reader would not be.
            var assay = Assay() with { Instrument = new CvParam("MS", "", "Orbitrap Fusion Lumos", "") };
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }).Results[0];

            Assert.That(row["comment[instrument]"], Is.EqualTo("NT=Orbitrap Fusion Lumos;AC=MS:1002732"));
        }

        [Test]
        public void AnUnresolvableInstrumentNameIsWrittenWithoutAnAccession()
        {
            var assay = Assay() with { Instrument = new CvParam("MS", "", "Some Prototype Instrument", "") };
            var row = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }).Results[0];

            Assert.That(row["comment[instrument]"], Is.EqualTo("NT=Some Prototype Instrument"));
            Assert.That(row["comment[instrument]"], Does.Not.Contain("AC="),
                "the spec permits a name alone; inventing an accession would be worse");
        }

        // ---------------- document shape ----------------

        [Test]
        public void RowsPadToTheWidestModificationCount()
        {
            var many = Assay("a.raw") with
            {
                VariableModifications = new[] { Mod("Oxidation on M", "M", "35"), Mod("Acetyl on K", "K", "1") }
            };
            var few = Assay("b.raw") with { VariableModifications = Array.Empty<Modification>() };

            var document = SdrfBuilder.Build(new[]
            {
                new SdrfRowInput(Sample("Sample 1"), many),
                new SdrfRowInput(Sample("Sample 2"), few)
            });

            Assert.That(document.Header.IndexesOf("comment[modification parameters]").Count, Is.EqualTo(3));
            Assert.That(document.Results[1].All("comment[modification parameters]").Skip(1),
                Is.All.EqualTo("not applicable"));
            Assert.That(document.Results.Select(r => r.Cells.Count).Distinct().Count(), Is.EqualTo(1),
                "a ragged document is what the reader had to be taught to tolerate; do not author one");
        }

        [Test]
        public void ColumnsAreInSpecificationBlockOrder()
        {
            var document = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) },
                new SdrfBuilderOptions
                {
                    ProteomeXchangeAccession = "PXD000070",
                    Software = new CvParam("MS", "MS:1002826", "MetaMorpheus", ""),
                    SoftwareVersion = "1.0.0"
                });
            var names = document.Header.ToList();

            Assert.That(names.IndexOf("source name"), Is.EqualTo(0));
            Assert.That(names.IndexOf("characteristics[organism]"),
                Is.LessThan(names.IndexOf("assay name")));
            Assert.That(names.IndexOf("assay name"), Is.LessThan(names.IndexOf("comment[instrument]")));
            Assert.That(SdrfValidator.Validate(document).Warnings.Any(w => w.Rule == "ColumnOrdering"),
                Is.False);
        }

        [Test]
        public void ProvenanceTiesOurAssayParametersToTheirSamples()
        {
            var document = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), Assay()) },
                new SdrfBuilderOptions
                {
                    ProteomeXchangeAccession = "PXD000070",
                    Software = new CvParam("MS", "MS:1002826", "MetaMorpheus", ""),
                    SoftwareVersion = "1.0.7"
                });
            var row = document.Results[0];

            Assert.That(row["comment[proteomexchange accession number]"], Is.EqualTo("PXD000070"));
            Assert.That(row["comment[software]"], Does.Contain("AC=MS:1002826").And.Contain("VV=1.0.7"));
            Assert.That(row["comment[sdrf version]"], Is.EqualTo("v1.1.0"));
        }

        [Test]
        public void TheBuiltDocumentRoundTrips()
        {
            var document = SdrfBuilder.Build(new[]
            {
                new SdrfRowInput(Sample("Sample 1"), Assay("a.raw")),
                new SdrfRowInput(Sample("Sample 2"), Assay("b.raw"))
            });

            string path = System.IO.Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"built_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(path);
                var reloaded = new SdrfDocument(path);

                Assert.That(reloaded.Header.ToList(), Is.EqualTo(document.Header.ToList()));
                Assert.That(reloaded.Results.Count, Is.EqualTo(2));
                Assert.That(reloaded.Results[0].Cells, Is.EqualTo(document.Results[0].Cells));
                Assert.That(SdrfValidator.Validate(reloaded).IsValid, Is.True);
            }
            finally
            {
                if (System.IO.File.Exists(path)) System.IO.File.Delete(path);
            }
        }

        /// <summary>
        /// A caller that supplies no characteristics at all still gets every required column.
        ///
        /// This is the one gap the rest of this fixture could not see, because Sample() has always
        /// populated organism part — so every other test here proved the builder writes a column it
        /// was handed, never that it writes one it was not. MetaMorpheus, which populates nothing,
        /// produced documents with the column ABSENT for as long as that was true.
        ///
        /// Absent is strictly worse than empty. SdrfCoverage reports the fill rate of columns that
        /// exist, so it cannot report a column that is not there; the document reads as complete
        /// precisely where it says least.
        /// </summary>
        [Test]
        public void RequiredCharacteristicsAreEmittedEvenWhenTheCallerSuppliesNone()
        {
            var sample = Sample() with { Characteristics = new Dictionary<string, CvParam>() };
            var options = new SdrfBuilderOptions { RequireSampleMetadata = false };

            var document = SdrfBuilder.Build(new[] { new SdrfRowInput(sample, Assay()) }, options);

            Assert.That(document.Header.Contains("characteristics[organism part]"), Is.True,
                "SDRF requires this column of every document; a caller that supplies no value must " +
                "still get the column, or nothing downstream can see that it says nothing.");
            Assert.That(document.Results[0]["characteristics[organism part]"], Is.EqualTo("not available"),
                "The specification permits the reserved word here, and a fifth of the curated corpus " +
                "uses one. Empty and honest beats absent.");
            Assert.That(SdrfValidator.Validate(document).IsValid, Is.True);
        }

        /// <summary>
        /// The strict path is unchanged: asking for sample metadata and not supplying it still
        /// throws rather than quietly writing a reserved word.
        /// </summary>
        [Test]
        public void RequiringSampleMetadataStillRefusesAnUnsuppliedRequiredCharacteristic()
        {
            var sample = Sample() with { Characteristics = new Dictionary<string, CvParam>() };

            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(sample, Assay()) }),
                Throws.TypeOf<MzLibException>());
        }

        [Test]
        public void EmptyInputIsRejected()
        {
            Assert.That(() => SdrfBuilder.Build(Array.Empty<SdrfRowInput>()),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(null), Throws.TypeOf<ArgumentNullException>());
        }

        /// <summary>
        /// RequireSampleMetadata = false exists so that a finished search is never thrown away, so
        /// nothing under it may throw. SdrfCell.ToCell rightly refuses a term with neither a name nor
        /// an accession, and three cell writers could hand it one: an instrument, a digestion agent
        /// and a modification that each name nothing.
        /// </summary>
        [Test]
        public void ALenientBuildDoesNotThrowOnATermThatNamesNothing()
        {
            var assay = Assay() with
            {
                Instrument = new CvParam("MS", "", "", ""),
                CleavageAgent = new Protease("", CleavageSpecificity.Full, "", "", new List<DigestionMotif>()),
                FixedModifications = Array.Empty<Modification>(),
                VariableModifications = new[] { new Modification() }
            };
            var options = new SdrfBuilderOptions { RequireSampleMetadata = false };

            SdrfDocument document = null;
            Assert.That(() => document = SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }, options),
                Throws.Nothing);

            Assert.That(document.Results[0]["comment[instrument]"], Is.EqualTo("not available"));
            Assert.That(document.Results[0]["comment[cleavage agent details]"], Is.EqualTo("not available"));
            Assert.That(document.Results[0]["comment[modification parameters]"], Is.EqualTo("not available"));
        }

        /// <summary>
        /// And the strict path still refuses it — but as the MzLibException this class contracts for,
        /// not as an ArgumentException escaping a formatting helper.
        /// </summary>
        [Test]
        public void TheStrictPathRefusesATermThatNamesNothingWithTheUsualException()
        {
            var assay = Assay() with { Instrument = new CvParam("MS", "", "", "") };

            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(Sample(), assay) }),
                Throws.TypeOf<MzLibException>());
        }

        /// <summary>
        /// Build null-checks its argument, so the members it then dereferences are checked too.
        /// Without this, a null Assay or Characteristics surfaced as a NullReferenceException thrown
        /// from inside a LINQ lambda, naming neither the row nor what was missing.
        /// </summary>
        [Test]
        public void ARowWithANullPartIsRejectedAsAnArgumentException()
        {
            var sample = Sample();
            var assay = Assay();

            Assert.That(() => SdrfBuilder.Build(new SdrfRowInput[] { null }),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(null, assay) }),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(sample, null) }),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(sample with { Characteristics = null }, assay) }),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(sample, assay with { FixedModifications = null }) }),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => SdrfBuilder.Build(new[] { new SdrfRowInput(sample, assay with { VariableModifications = null }) }),
                Throws.TypeOf<ArgumentException>());
        }
    }
}
