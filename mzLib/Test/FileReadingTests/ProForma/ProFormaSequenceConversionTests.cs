using System.Linq;
using NUnit.Framework;
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
    }
}
