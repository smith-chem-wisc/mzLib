using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Reading an isobaric kit from a record's text (sdrf PXR-F15 P1). Measured first on the 41 curated isobaric deposits
    /// with a cached PRIDE record (sdrf results/isobaric-p0-2026-11-16): a record names the kit with its size in 17 and is
    /// wrong in 1; PRIDE's quantification method only ever names the family; and a record whose only kit word is
    /// <c>TMT6plex</c> was a 10-plex in 3 of 12, because TMT6plex is also Unimod's name for the 6, 10 and 11-plex reagent.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfIsobaricKit
    {
        [TestCase("Peptides were labelled with TMT10plex reagents (Thermo).", IsobaricMassTagType.TMT10)]
        [TestCase("labeled using a TMT 10-plex kit", IsobaricMassTagType.TMT10)]
        [TestCase("a 10-plex TMT experiment", IsobaricMassTagType.TMT10)]
        [TestCase("TMT11plex Isobaric Label Reagent Set", IsobaricMassTagType.TMT11)]
        [TestCase("TMT 11-plex", IsobaricMassTagType.TMT11)]
        [TestCase("TMTpro 16plex Label Reagent Set", IsobaricMassTagType.TMT16)]
        [TestCase("labelled with TMTpro-16", IsobaricMassTagType.TMT16)]
        [TestCase("TMTpro 18plex", IsobaricMassTagType.TMT18)]
        [TestCase("TMTpro18 reagents", IsobaricMassTagType.TMT18)]
        [TestCase("iTRAQ 8-plex", IsobaricMassTagType.iTRAQ8)]
        [TestCase("labeled with an iTRAQ Reagent-8plex Multiplex Kit (AB Sciex)", IsobaricMassTagType.iTRAQ8)]
        [TestCase("iTRAQ8plex", IsobaricMassTagType.iTRAQ8)]
        [TestCase("4plex iTRAQ", IsobaricMassTagType.iTRAQ4)]
        [TestCase("iTRAQ 4-plex reagents", IsobaricMassTagType.iTRAQ4)]
        [TestCase("TMT10plex labelling; TMT6plex (K) and TMT6plex (N-term) were set as fixed modifications", IsobaricMassTagType.TMT10,
            Description = "TMT6plex beside a sized kit is the modification's name, not a second kit")]
        [TestCase("TMT reagent 10 plex; mobile phase A (95% water, 5% dimethylsulfoxide (DMSO))", IsobaricMassTagType.TMT10,
            Description = "PXD017201: a solvent, not dimethyl labelling")]
        [TestCase("the pellet resuspended in label-free buffer and labelled with TMT10plex", IsobaricMassTagType.TMT10,
            Description = "PXD013753: a buffer's name, not label-free quantification")]
        public void ARecordThatNamesTheKitWithItsSizeGivesThatKit(string text, IsobaricMassTagType kit)
        {
            var reading = SdrfIsobaricKit.Read(text);

            Assert.That(reading.Kit, Is.EqualTo(kit));
            Assert.That(reading.Labelled, Is.True);
            Assert.That(text, Does.Contain(reading.Quote), "the quote is the record's own words");
        }

        [TestCase("Fixed modifications: TMT6plex (K), TMT6plex (N-term)", Description = "Unimod's name for the 6, 10 and 11-plex reagent")]
        [TestCase("labelled with TMT 6-plex reagents", Description = "measured: a record whose only kit word is TMT6plex was a 10-plex in 3 of 12")]
        [TestCase("TMTpro labelling", Description = "16 or 18")]
        [TestCase("tandem mass tag labelling")]
        [TestCase("iTRAQ labelling")]
        [TestCase("TMT10plex for the first set and TMT11plex for the second", Description = "two kits named: neither is the deposit's")]
        [TestCase("labeled with the corresponding TMT 10-plex/ 6-plex reagent", Description = "PXD009243: a 6-plex written apart is a kit beside the 10-plex")]
        [TestCase("to compare label-free quantification with TMT10plex labelling", Description = "PXD004682's shape: some runs are label free")]
        [TestCase("SILAC labelling, and TMT 10plex labelling kits for the TMT runs", Description = "PXD014565: label-free, SILAC and TMT runs in one deposit")]
        [TestCase("18O labelling of the native samples, and a TMT10plex kit", Description = "PXD030245")]
        public void ARecordThatNamesOnlyTheFamilyIsLabelledButNamesNoKit(string text)
        {
            var reading = SdrfIsobaricKit.Read(text);

            Assert.That(reading.Kit, Is.Null);
            Assert.That(reading.Labelled, Is.True);
        }

        [TestCase(null)]
        [TestCase("")]
        [TestCase("Label-free quantification with MaxLFQ.")]
        [TestCase("Samples were measured at 10 plex-like intervals", Description = "no TMT or iTRAQ word")]
        [TestCase("the mtmt10 gene", Description = "a word that only contains the letters")]
        public void TextNamingNoIsobaricLabelIsNotLabelled(string? text)
        {
            var reading = SdrfIsobaricKit.Read(text);

            Assert.That(reading.Labelled, Is.False);
            Assert.That(reading.Kit, Is.Null);
        }
    }
}
