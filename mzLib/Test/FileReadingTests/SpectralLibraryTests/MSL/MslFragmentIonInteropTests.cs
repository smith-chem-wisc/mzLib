using System.Collections.Generic;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Guards the mapping between the spectral library's own fragment model (<see cref="MslFragmentIon"/>)
/// and the search/scoring domain type (<see cref="MatchedFragmentIon"/>). The two are deliberately
/// separate types bridged only by <see cref="MslLibraryEntry.ToLibrarySpectrum"/> /
/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/>; these tests round-trip a fragment through that seam
/// and assert every mapped field survives, so the two models cannot silently drift apart. They also pin
/// the one field (<see cref="MslFragmentIon.ExcludeFromQuant"/>) that intentionally does NOT survive the
/// domain projection.
/// </summary>
[TestFixture]
public class MslFragmentIonInteropTests
{
    private static MslLibraryEntry MakeEntry(IEnumerable<MslFragmentIon> fragments) => new()
    {
        FullSequence = "PEPTIDEK",
        BaseSequence = "PEPTIDEK",
        PrecursorMz = 456.789,
        ChargeState = 2,
        RetentionTime = 33.3,
        IsDecoy = false,
        MatchedFragmentIons = new List<MslFragmentIon>(fragments)
    };

    /// <summary>
    /// A deliberately diverse fragment set — terminal ion, terminal ion with neutral loss at charge 2,
    /// internal ion (bIy[3-6]), and a diagnostic ion — round-trips MslFragmentIon -> MatchedFragmentIon
    /// (via LibrarySpectrum) -> MslFragmentIon with every mapped field preserved. If a new mapped field is
    /// added to one side of the conversion but not the other, this test fails.
    /// </summary>
    [Test]
    public void RoundTripThroughMatchedFragmentIon_PreservesEveryMappedField()
    {
        var original = MakeEntry(new[]
        {
            // terminal b5, charge 1, no loss
            new MslFragmentIon { Mz = 100.1f, Intensity = 0.5f, ProductType = ProductType.b, FragmentNumber = 5, ResiduePosition = 5, Charge = 1 },
            // terminal y3 with a water loss, charge 2
            new MslFragmentIon { Mz = 200.2f, Intensity = 1.0f, ProductType = ProductType.y, FragmentNumber = 3, ResiduePosition = 3, Charge = 2, NeutralLoss = 18.0106 },
            // internal bIy spanning residues 3..6
            new MslFragmentIon { Mz = 300.3f, Intensity = 0.25f, ProductType = ProductType.b, SecondaryProductType = ProductType.y, FragmentNumber = 3, SecondaryFragmentNumber = 6, ResiduePosition = 3, Charge = 1 },
            // diagnostic ion
            new MslFragmentIon { Mz = 138.06f, Intensity = 0.1f, ProductType = ProductType.D, FragmentNumber = 138, ResiduePosition = 0, Charge = 1 },
        });

        // MslFragmentIon -> MatchedFragmentIon (LibrarySpectrum) -> MslFragmentIon
        LibrarySpectrum asDomain = original.ToLibrarySpectrum();
        MslLibraryEntry roundTripped = MslLibraryEntry.FromLibrarySpectrum(asDomain);

        Assert.That(roundTripped, Is.Not.Null);
        Assert.That(roundTripped.MatchedFragmentIons, Has.Count.EqualTo(original.MatchedFragmentIons.Count));

        for (int i = 0; i < original.MatchedFragmentIons.Count; i++)
        {
            MslFragmentIon before = original.MatchedFragmentIons[i];
            MslFragmentIon after = roundTripped.MatchedFragmentIons[i];

            // float -> double -> float is exactly lossless, so these compare exactly.
            Assert.Multiple(() =>
            {
                Assert.That(after.Mz, Is.EqualTo(before.Mz), $"Mz[{i}]");
                Assert.That(after.Intensity, Is.EqualTo(before.Intensity), $"Intensity[{i}]");
                Assert.That(after.Charge, Is.EqualTo(before.Charge), $"Charge[{i}]");
                Assert.That(after.ProductType, Is.EqualTo(before.ProductType), $"ProductType[{i}]");
                Assert.That(after.SecondaryProductType, Is.EqualTo(before.SecondaryProductType), $"SecondaryProductType[{i}]");
                Assert.That(after.FragmentNumber, Is.EqualTo(before.FragmentNumber), $"FragmentNumber[{i}]");
                Assert.That(after.SecondaryFragmentNumber, Is.EqualTo(before.SecondaryFragmentNumber), $"SecondaryFragmentNumber[{i}]");
                Assert.That(after.ResiduePosition, Is.EqualTo(before.ResiduePosition), $"ResiduePosition[{i}]");
                Assert.That(after.NeutralLoss, Is.EqualTo(before.NeutralLoss), $"NeutralLoss[{i}]");
                // Derived flags must agree as a consequence of the above mapping.
                Assert.That(after.IsInternalFragment, Is.EqualTo(before.IsInternalFragment), $"IsInternalFragment[{i}]");
                Assert.That(after.IsDiagnosticIon, Is.EqualTo(before.IsDiagnosticIon), $"IsDiagnosticIon[{i}]");
            });
        }
    }

    /// <summary>
    /// Documents the one deliberate gap in the mapping: ExcludeFromQuant is an MSL-only concept with no
    /// home in MatchedFragmentIon / LibrarySpectrum, so it is dropped by the domain projection. (It DOES
    /// round-trip through the binary .msl file — see TestMslCoverage.) Pinning this prevents someone from
    /// "fixing" the perceived data loss by smuggling the flag through the domain type.
    /// </summary>
    [Test]
    public void ExcludeFromQuant_IsNotPreservedThroughLibrarySpectrumProjection()
    {
        var entry = MakeEntry(new[]
        {
            new MslFragmentIon { Mz = 100.1f, Intensity = 0.5f, ProductType = ProductType.b, FragmentNumber = 5, ResiduePosition = 5, Charge = 1, ExcludeFromQuant = true },
        });

        MslLibraryEntry roundTripped = MslLibraryEntry.FromLibrarySpectrum(entry.ToLibrarySpectrum());

        Assert.That(entry.MatchedFragmentIons[0].ExcludeFromQuant, Is.True, "set true on the original library fragment");
        Assert.That(roundTripped.MatchedFragmentIons[0].ExcludeFromQuant, Is.False,
            "ExcludeFromQuant has no representation in MatchedFragmentIon/LibrarySpectrum, so the domain projection intentionally drops it; it round-trips only through the binary .msl format.");
    }
}
