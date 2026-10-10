using System.Collections.Generic;
using System.Collections.Immutable;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using Readers.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// ProForma as a registered format of the shared <see cref="SequenceConversionService"/>
    /// (Omics/SequenceConversion). Covers the flat-subset boundary — the point at which a lossless
    /// Layer-1 term stops fitting the flat <see cref="CanonicalSequence"/> model — and the
    /// cross-format conversions that registration buys.
    /// </summary>
    [TestFixture]
    internal class ProFormaSequenceConversionTests
    {
        private SequenceConversionService _service = null!;

        [SetUp]
        public void SetUp()
        {
            // A private service rather than Default, so registration order cannot leak between fixtures.
            _service = new SequenceConversionService();
            _service.RegisterParser(MzLibSequenceParser.Instance);
            _service.RegisterSerializer(MzLibSequenceSerializer.Instance);
            ProFormaSequenceConversion.RegisterWith(_service);
        }

        [Test]
        public void Register_MakesProFormaASourceAndTargetFormat()
        {
            Assert.That(_service.AvailableSourceFormats, Does.Contain("ProForma"));
            Assert.That(_service.AvailableTargetFormats, Does.Contain("ProForma"));
        }

        [Test]
        public void Schema_DescribesProFormaBracketAndTerminalSyntax()
        {
            var schema = ProFormaSequenceParser.Instance.Schema;

            Assert.That(schema.FormatName, Is.EqualTo("ProForma"));
            Assert.That(schema, Is.SameAs(ProFormaSequenceFormatSchema.Instance));
            Assert.That(schema.ModOpenBracket, Is.EqualTo('['));
            Assert.That(schema.ModCloseBracket, Is.EqualTo(']'));
            Assert.That(schema.SupportsNTerminalMods, Is.True);
            Assert.That(schema.SupportsCTerminalMods, Is.True);
        }

        [Test]
        public void Parse_AccessionTag_PopulatesCanonicalSequence()
        {
            var canonical = _service.Parse("PEPC[UNIMOD:4]TIDE", "ProForma");

            Assert.That(canonical, Is.Not.Null);
            Assert.That(canonical!.Value.BaseSequence, Is.EqualTo("PEPCTIDE"));
            Assert.That(canonical.Value.ModificationCount, Is.EqualTo(1));

            var mod = canonical.Value.GetModificationAt(3);
            Assert.That(mod, Is.Not.Null);
            Assert.That(mod!.Value.UnimodId, Is.EqualTo(4));
            Assert.That(mod.Value.TargetResidue, Is.EqualTo('C'));
            Assert.That(mod.Value.OriginalRepresentation, Is.EqualTo("UNIMOD:4"));
        }

        [Test]
        public void Parse_TerminalDescriptors_MapToTerminalPositions()
        {
            var canonical = _service.Parse("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]", "ProForma");

            Assert.That(canonical, Is.Not.Null);
            Assert.That(canonical!.Value.NTerminalModification?.UnimodId, Is.EqualTo(1));
            Assert.That(canonical.Value.CTerminalModification?.UnimodId, Is.EqualTo(2));
            Assert.That(canonical.Value.ResidueModifications.Count(), Is.Zero);
        }

        [TestCase("PEPTIDE-[UNIMOD:2][UNIMOD:35]", "Multiple C-terminal modifications (2 found).")]
        [TestCase("[UNIMOD:1][UNIMOD:35]-PEPTIDE", "Multiple N-terminal modifications (2 found).")]
        public void Parse_StackedTerminalModifications_AreIncompatibleNotLost(string input, string item)
        {
            var warnings = new ConversionWarnings();
            var canonical = ProFormaSequenceParser.Instance.Parse(input, warnings, SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
            Assert.That(warnings.IncompatibleItems, Does.Contain(item));
            Assert.That(warnings.Errors, Has.Some.Contains("ProForma 2.1"));
            Assert.That(() => ProFormaSequenceParser.Instance.Parse(input), Throws.TypeOf<SequenceConversionException>());
            Assert.That(ProFormaSequenceParser.Instance.CanParse(input), Is.False);
        }

        [Test]
        public void Parse_ResidueAfterCTerminalModification_IsInvalidSequence()
        {
            var warnings = new ConversionWarnings();
            var canonical = ProFormaSequenceParser.Instance.Parse("PEPTIDE-[UNIMOD:2]K", warnings, SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
            Assert.That(warnings.Errors, Has.Some.Contains("Unexpected content at position 18 after the C-terminal modification."));
        }

        [Test]
        public void RoundTrip_ProFormaToProForma_PreservesAccessionNotation()
        {
            var proForma = _service.Convert("PEPC[UNIMOD:4]TIDE", "ProForma", "ProForma");

            Assert.That(proForma, Is.EqualTo("PEPC[UNIMOD:4]TIDE"));
        }

        [Test]
        public void Convert_MzLibToProForma_EmitsStandardAccession()
        {
            // The headline interop step: a tool-specific mzLib name becomes a community accession,
            // via Nic's mzLib parser + our ProForma serializer, with no bespoke converter.
            var proForma = _service.Convert("PEPC[Common Fixed:Carbamidomethyl on C]TIDE", "mzLib", "ProForma");

            Assert.That(proForma, Is.EqualTo("PEPC[UNIMOD:4]TIDE"));
        }

        [Test]
        public void Parse_MassTag_KeepsTheMassShift()
        {
            var canonical = _service.Parse("PEPM[+15.9949]TIDE", "ProForma");

            Assert.That(canonical, Is.Not.Null);
            var mod = canonical!.Value.GetModificationAt(3);
            Assert.That(mod, Is.Not.Null);
            Assert.That(mod!.Value.MonoisotopicMass, Is.EqualTo(15.9949).Within(1e-6));
        }

        // The flat-subset boundary. Each of these is valid ProForma 2.0 that Layer 1 round-trips
        // losslessly, but which CanonicalSequence cannot represent - so it must fail loud, not
        // silently drop data.
        [TestCase("PEPTIDEC[XLMOD:02001#XL1]KAC[#XL1]R", "tag groups", TestName = "TagGroup_IsUnrepresentable")]
        [TestCase("{Glycan:Hex}PEPTIDE", "labile modifications", TestName = "Labile_IsUnrepresentable")]
        [TestCase("[UNIMOD:21]?PEPTIDE", "unlocalized modifications", TestName = "Unlocalized_IsUnrepresentable")]
        [TestCase("<[UNIMOD:4]@C>PEPCTIDE", "global modifications", TestName = "Global_IsUnrepresentable")]
        [TestCase("PROT(EOSFORMS)[+19.0523]ISK", "position ranges", TestName = "Range_IsUnrepresentable")]
        // A single-residue ambiguous region: the range check (start != end) does not fire, so this is
        // the case that actually reaches the ambiguity guard. A multi-residue "(?DQ)" is caught too,
        // but reports as "position ranges" - whichever guard sees it first.
        [TestCase("PROT(?D)NGTWEM[UNIMOD:35]ESNENFEGYMK", "sequence ambiguity", TestName = "Ambiguity_IsUnrepresentable")]
        public void Parse_ConstructBeyondTheFlatModel_FailsLoud(string proForma, string expectedReason)
        {
            // Layer 1 handles it losslessly...
            Assert.That(() => ProFormaReader.Read(proForma), Throws.Nothing);

            // ...but the canonical model must refuse it rather than lose it.
            var warnings = new ConversionWarnings();
            Assert.That(() => _service.Parse(proForma, "ProForma", warnings),
                Throws.TypeOf<SequenceConversionException>());

            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
            Assert.That(warnings.IncompatibleItems, Has.Some.Contains(expectedReason));
        }

        [Test]
        public void Parse_ConstructBeyondTheFlatModel_ReturnsNullWhenAsked()
        {
            var warnings = new ConversionWarnings();
            var canonical = _service.Parse("{Glycan:Hex}PEPTIDE", "ProForma", warnings,
                SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.HasFatalError, Is.True);
        }

        [Test]
        public void Parse_InvalidProForma_IsReportedAsInvalidSequence()
        {
            var warnings = new ConversionWarnings();
            var canonical = _service.Parse("PEP[UNIMOD:4TIDE", "ProForma", warnings,
                SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
        }

        [Test]
        public void MultipleDescriptors_KeepTheAccessionAndWarn()
        {
            var warnings = new ConversionWarnings();
            var canonical = _service.Parse("PEPC[Carbamidomethyl|UNIMOD:4]TIDE", "ProForma", warnings);

            Assert.That(canonical, Is.Not.Null);
            var mod = canonical!.Value.GetModificationAt(3);
            Assert.That(mod!.Value.UnimodId, Is.EqualTo(4), "the accession is the most resolvable descriptor");
            Assert.That(warnings.HasWarnings, Is.True, "dropping a descriptor is lossy and must be said out loud");
        }

        [Test]
        public void CanParse_LeavesMzLibStringsToTheMzLibParser()
        {
            // Auto-detect must not let ProForma steal the mzLib format's own syntax.
            Assert.That(ProFormaSequenceParser.Instance.CanParse("PEPC[Common Fixed:Carbamidomethyl on C]TIDE"), Is.False);
            Assert.That(ProFormaSequenceParser.Instance.CanParse("PEPTIDE"), Is.False);

            Assert.That(ProFormaSequenceParser.Instance.CanParse("PEPC[UNIMOD:4]TIDE"), Is.True);
            Assert.That(ProFormaSequenceParser.Instance.CanParse("PEPM[+15.9949]TIDE"), Is.True);
        }

        [Test]
        public void CanParse_RejectsConstructsItCannotRepresent()
        {
            // Valid ProForma, but not flat - claiming it in auto-detect would guarantee a failed parse.
            Assert.That(ProFormaSequenceParser.Instance.CanParse("{Glycan:Hex}PEPTIDE"), Is.False);
        }

        [Test]
        public void RegisterWithDefault_IsIdempotent()
        {
            ProFormaSequenceConversion.RegisterWithDefault();
            ProFormaSequenceConversion.RegisterWithDefault(); // second call must be a no-op, not a double registration

            Assert.That(SequenceConversionService.Default.AvailableSourceFormats, Does.Contain("ProForma"));
            Assert.That(SequenceConversionService.Default.AvailableTargetFormats, Does.Contain("ProForma"));
        }

        [Test]
        public void Serialize_TerminalMods_RoundTripThroughTheService()
        {
            var proForma = _service.Convert("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]", "ProForma", "ProForma");

            Assert.That(proForma, Is.EqualTo("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]"));
        }

        [Test]
        public void Serialize_UnresolvableName_FallsBackToANameDescriptor()
        {
            // The lookup cannot resolve this, so the serializer must keep the name rather than
            // invent an accession or drop the modification.
            var proForma = _service.Convert("PEPC[NotARealModification]TIDE", "ProForma", "ProForma");

            Assert.That(proForma, Is.EqualTo("PEPC[NotARealModification]TIDE"));
        }

        [Test]
        public void WithLookup_UsesTheSuppliedLookup()
        {
            // The default lookup resolves "Carbamidomethyl" to an mzLib mod and emits [UNIMOD:4]
            // (see Convert_MzLibToProForma). A lookup that resolves nothing must leave the name alone.
            var serializer = ProFormaSequenceSerializer.WithLookup(NeverResolvesLookup.Instance);
            var canonical = _service.Parse("PEPC[Carbamidomethyl]TIDE", "ProForma");

            Assert.That(serializer.ModificationLookup, Is.SameAs(NeverResolvesLookup.Instance));
            Assert.That(serializer.Serialize(canonical!.Value), Is.EqualTo("PEPC[Carbamidomethyl]TIDE"));
        }

        [Test]
        public void Serialize_MassOnlyMod_FallsBackToAMassDescriptor()
        {
            // No name, no accession, but a mass the lookup cannot match: ProForma can still say the mass.
            var sequence = SequenceWithSingleMod(CanonicalModification.AtResidue(
                residueIndex: 3, targetResidue: 'C', originalRepresentation: string.Empty, mass: 999.999));

            var proForma = ProFormaSequenceSerializer.WithLookup(NeverResolvesLookup.Instance).Serialize(sequence);

            Assert.That(proForma, Is.EqualTo("PEPC[+999.999]TIDE"));
        }

        [Test]
        public void Serialize_ModWithNoDescriptor_FailsLoud()
        {
            // Nothing to say about this modification at all - dropping it silently would lose data.
            var sequence = SequenceWithSingleMod(CanonicalModification.AtResidue(
                residueIndex: 3, targetResidue: 'C', originalRepresentation: string.Empty));

            var warnings = new ConversionWarnings();
            Assert.That(() => ProFormaSequenceSerializer.WithLookup(NeverResolvesLookup.Instance)
                    .Serialize(sequence, warnings),
                Throws.TypeOf<SequenceConversionException>());

            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
            Assert.That(warnings.IncompatibleItems, Has.Some.Contains("cannot be written as ProForma"));
        }

        [Test]
        public void Serialize_ModWithNoDescriptor_IsDroppedWhenAsked()
        {
            var sequence = SequenceWithSingleMod(CanonicalModification.AtResidue(
                residueIndex: 3, targetResidue: 'C', originalRepresentation: string.Empty));

            var warnings = new ConversionWarnings();
            var proForma = ProFormaSequenceSerializer.WithLookup(NeverResolvesLookup.Instance)
                .Serialize(sequence, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

            Assert.That(proForma, Is.EqualTo("PEPCTIDE"), "the mod is dropped, the sequence survives");
            Assert.That(warnings.HasWarnings, Is.True, "a silent drop is exactly what we must not do");
        }

        [Test]
        public void Serialize_EmptyBaseSequence_IsInvalid()
        {
            var warnings = new ConversionWarnings();
            var proForma = ProFormaSequenceSerializer.Instance.Serialize(
                CanonicalSequence.Empty, warnings, SequenceConversionHandlingMode.ReturnNull);

            Assert.That(proForma, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
        }

        [Test]
        public void Parse_NullOrEmpty_IsInvalidSequence()
        {
            // Straight at the parser: the service short-circuits empty input before the parser sees it.
            var warnings = new ConversionWarnings();
            var canonical = ProFormaSequenceParser.Instance.Parse(
                "   ", warnings, SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.InvalidSequence));
        }

        [Test]
        public void CanParse_ProFormaLookalikeThatDoesNotParse_IsFalse()
        {
            // Carries a ProForma token but is malformed - must not be claimed by auto-detection.
            Assert.That(ProFormaSequenceParser.Instance.CanParse("PEPC[UNIMOD:4TIDE"), Is.False);
        }

        [Test]
        public void Parse_TwoModsOnOneResidue_IsRejected()
        {
            // Valid ProForma; the flat canonical model allows only one modification per position.
            var warnings = new ConversionWarnings();
            var canonical = _service.Parse("PEPC[UNIMOD:4][UNIMOD:35]TIDE", "ProForma", warnings,
                SequenceConversionHandlingMode.ReturnNull);

            Assert.That(canonical, Is.Null);
            Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
            Assert.That(warnings.IncompatibleItems, Has.Some.Contains("Multiple modifications at residue index 3"));
        }

        [Test]
        public void Parse_NonUnimodOntology_HasNoUnimodId()
        {
            var canonical = _service.Parse("PEPC[MOD:00397]TIDE", "ProForma");

            var mod = canonical!.Value.GetModificationAt(3);
            Assert.That(mod!.Value.OriginalRepresentation, Is.EqualTo("MOD:00397"));
            Assert.That(mod.Value.UnimodId, Is.Null, "only UNIMOD accessions populate UnimodId");
        }

        [Test]
        public void Instance_DefaultLookup_IsMetaMorpheusFirst()
        {
            // The MetaMorpheus-only lookup never resolved UniProt mods (#1401); Global alone moved mass-only mods
            // to a different isobaric Unimod record. MetaMorpheus first, then Global, does neither.
            Assert.That(ProFormaSequenceSerializer.Instance.ModificationLookup,
                Is.SameAs(ProFormaSequenceSerializer.MetaMorpheusFirstLookup.Instance));
        }

        [TestCase("PEN[+0.984016]K", "MassShift", "PEN[UNIMOD:7]K")]
        [TestCase("PEQ[+0.984016]K", "MassShift", "PEQ[UNIMOD:7]K")]
        [TestCase("PEF[+15.9949]K", "ProForma", "PEF[UNIMOD:35]K")]
        [TestCase("PEK[+42.0470]K", "ProForma", "PEK[UNIMOD:37]K")]
        [TestCase("PEK[+57.0215]K", "ProForma", "PEK[UNIMOD:4]K")]
        public void MassOnlyMod_KeepsItsMetaMorpheusIdentity(string sequence, string sourceFormat, string expected)
        {
            // With Global alone these became Asn->Asp, Gln->Glu, Phe->Tyr, Propyl and Gly: same mass, other identity.
            Assert.That(SequenceConversionService.Default.Convert(sequence, sourceFormat, "ProForma"), Is.EqualTo(expected));
        }

        /// <summary>
        /// Every modification mzLib ships as a MetaMorpheus modification, placed where its location
        /// restriction allows, as an mzLib full sequence.
        /// </summary>
        private static IEnumerable<TestCaseData> MetaMorpheusModFullSequences()
        {
            foreach (var mod in Mods.MetaMorpheusModifications)
            {
                char residue = (mod.Target?.ToString() ?? "X").FirstOrDefault(char.IsUpper);
                if (residue == default || residue == 'X') residue = 'A';
                string tag = $"[{mod.ModificationType}:{mod.IdWithMotif}]";
                string restriction = mod.LocationRestriction ?? "";
                string fullSequence = restriction.Contains("N-terminal") || restriction.Contains("5'-terminal") ? $"{tag}{residue}EPK"
                    : restriction.Contains("C-terminal") || restriction.Contains("3'-terminal") ? $"PEP{residue}-{tag}"
                    : $"PE{residue}{tag}K";
                yield return new TestCaseData(fullSequence).SetArgDisplayNames($"{mod.ModificationType}:{mod.IdWithMotif}");
            }
        }

        [TestCaseSource(nameof(MetaMorpheusModFullSequences))]
        public void MzLibToProForma_MetaMorpheusMod_IsUnchangedByTheDefault(string fullSequence)
        {
            // The default moved from MzLibModificationLookup to MetaMorpheusFirstLookup. A MetaMorpheus
            // modification must still resolve to the entry it resolved to before, so its ProForma is unchanged.
            var before = new SequenceConverter(MzLibSequenceParser.Instance,
                ProFormaSequenceSerializer.WithLookup(MzLibModificationLookup.Instance));
            var after = new SequenceConverter(MzLibSequenceParser.Instance, ProFormaSequenceSerializer.Instance);

            Assert.That(after.Convert(fullSequence, mode: SequenceConversionHandlingMode.ReturnNull),
                Is.EqualTo(before.Convert(fullSequence, mode: SequenceConversionHandlingMode.ReturnNull)));
        }

        [TestCase("[UniProt:N-acetylserine on S]SEQK", "[UNIMOD:1]-SEQK")]
        [TestCase("PEPK[UniProt:N6,N6-dimethyllysine on K]R", "PEPK[UNIMOD:36]R")]
        [TestCase("PEPM[Common Variable:Oxidation on M]K", "PEPM[UNIMOD:35]K")]
        public void MzLibToProForma_UniProtMod_IsWrittenAsItsUnimodAccession(string fullSequence, string expected)
        {
            Assert.That(_service.Convert(fullSequence, "mzLib", "ProForma"), Is.EqualTo(expected));
        }

        private static CanonicalSequence SequenceWithSingleMod(CanonicalModification mod) =>
            new("PEPCTIDE", ImmutableArray.Create(mod), "test");

        /// <summary>A lookup that resolves nothing, so serializer fallbacks are reachable in a test.</summary>
        private sealed class NeverResolvesLookup : IModificationLookup
        {
            public static NeverResolvesLookup Instance { get; } = new();

            public string Name => "never-resolves";

            public CanonicalModification? TryResolve(CanonicalModification mod) => null;

            public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null,
                ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null) => null;
        }
    }
}
