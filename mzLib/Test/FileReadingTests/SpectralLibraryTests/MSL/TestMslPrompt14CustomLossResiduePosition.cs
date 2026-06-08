using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.IO;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 tests documenting and verifying the custom-neutral-loss ResiduePosition
/// limitation described in the Prompt 14 code review.
///
/// These tests:
///   1. Document the known limitation precisely (ResiduePosition = 0 after round-trip
///      for custom-loss fragments) with clear, explanatory failure messages.
///   2. Confirm the neutral-loss mass itself IS correctly recovered.
///   3. Confirm FragmentNumber is preserved and usable as a positional proxy.
///   4. Confirm named-loss fragments are NOT affected (parameterized over H2O, NH3, H3PO4, None).
///   5. Confirm multiple distinct custom losses in one entry all lose ResiduePosition.
///   6. Confirm a mixed entry (custom + named losses) only loses position for the custom fragment.
///   7. Confirm writing with ResiduePosition=0 is a clean no-op.
/// </summary>
[TestFixture]
public sealed class TestMslPrompt14CustomLossResiduePosition
{
    // ── Fixture paths ─────────────────────────────────────────────────────────

    private static readonly string OutputDir =
        Path.Combine(Path.GetTempPath(), "MslPrompt14Tests");

    [OneTimeSetUp]
    public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

    [OneTimeTearDown]
    public void OneTimeTearDown()
    {
        if (Directory.Exists(OutputDir))
            Directory.Delete(OutputDir, recursive: true);
    }

    private static string TempPath(string name) =>
        Path.Combine(OutputDir, name + ".msl");

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// <summary>
    /// Builds a minimal single-fragment library entry for round-trip tests.
    /// </summary>
    private static MslLibraryEntry EntryWithFragment(
        double neutralLoss,
        int residuePosition,
        ProductType productType = ProductType.y,
        int fragmentNumber = 3)
    {
        return new MslLibraryEntry
        {
            FullSequence = "PEPTIDE",
            BaseSequence = "PEPTIDE",
            PrecursorMz = 449.74,
            ChargeState = 2,
            MoleculeType = MslFormat.MoleculeType.Peptide,
            DissociationType = DissociationType.HCD,
            MatchedFragmentIons = new List<MslFragmentIon>
            {
                new MslFragmentIon
                {
                    Mz              = 300.0f,
                    Intensity       = 1.0f,
                    ProductType     = productType,
                    FragmentNumber  = fragmentNumber,
                    ResiduePosition = residuePosition,
                    Charge          = 1,
                    NeutralLoss     = neutralLoss
                }
            }
        };
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Core limitation: custom losses always yield ResiduePosition = 0
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Documents the known binary-format limitation: a custom-loss fragment written with
    /// ResiduePosition = 7 must be read back with ResiduePosition = 0, because the
    /// MslFragmentRecord.ResiduePosition field is repurposed as ExtAnnotationIdx at write
    /// time, and the original value is irrecoverably lost.
    ///
    /// This test is intentionally named "Documents" rather than asserting correct behaviour
    /// to signal to future maintainers that this is a known limitation, not an accident.
    /// See MslFragmentIon.ResiduePosition XML doc and the Prompt 14 review for the
    /// format-v4 upgrade path.
    /// </summary>
    [Test]
    public void CustomLoss_Documents_ResiduePositionLostAfterRoundTrip()
    {
        const double customLoss = -203.0794; // HexNAc — not a named loss
        const int originalResidPos = 7;
        string path = TempPath("residue_lost_custom");

        MslWriter.Write(path, new[] { EntryWithFragment(customLoss, originalResidPos) });

        using MslLibraryData lib = MslReader.Load(path);
        int recovered = lib.Entries[0].MatchedFragmentIons[0].ResiduePosition;

        Assert.That(recovered, Is.EqualTo(0),
            $"KNOWN LIMITATION (format v1–v3): The MslFragmentRecord.ResiduePosition " +
            $"field is repurposed as ExtAnnotationIdx for custom-loss fragments. The " +
            $"original value ({originalResidPos}) is irrecoverably lost in the binary " +
            $"record. ResiduePosition is always 0 after round-trip for custom-loss " +
            $"fragments. See MslFragmentIon.ResiduePosition XML doc and Prompt 14 " +
            $"review for the format-v4 upgrade path.");
    }

    /// <summary>
    /// Although ResiduePosition is lost, the neutral-loss mass itself must be correctly
    /// recovered from the extended annotation table. This confirms the annotation-index
    /// mechanism works for its intended purpose even when position data is sacrificed.
    /// </summary>
    [Test]
    public void CustomLoss_NeutralLossMass_IsCorrectlyRecovered()
    {
        const double customLoss = -203.0794;
        string path = TempPath("neutral_loss_recovered");

        MslWriter.Write(path, new[] { EntryWithFragment(customLoss, residuePosition: 5) });

        using MslLibraryData lib = MslReader.Load(path);
        double recovered = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

        Assert.That(recovered, Is.EqualTo(customLoss).Within(1e-4),
            "Custom neutral-loss mass must be correctly recovered even though " +
            "ResiduePosition is lost for the same fragment.");
    }

    /// <summary>
    /// FragmentNumber is preserved for custom-loss fragments and can serve as a positional
    /// proxy for consumers that need a cleavage-site index when ResiduePosition is
    /// unavailable (e.g. peak annotation labels, PTM site localization).
    /// </summary>
    [Test]
    public void CustomLoss_FragmentNumber_IsPreservedAsPositionalProxy()
    {
        const double customLoss = -203.0794;
        const int fragmentNumber = 4;
        string path = TempPath("fragment_number_proxy");

        MslWriter.Write(path, new[]
        {
            EntryWithFragment(customLoss, residuePosition: 7, fragmentNumber: fragmentNumber)
        });

        using MslLibraryData lib = MslReader.Load(path);
        int recovered = lib.Entries[0].MatchedFragmentIons[0].FragmentNumber;

        Assert.That(recovered, Is.EqualTo(fragmentNumber),
            "FragmentNumber must survive round-trip for custom-loss fragments and " +
            "can be used as a positional proxy when ResiduePosition is unavailable.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Named losses are NOT affected — parameterized
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// For named neutral losses (H2O, NH3, H3PO4, and no-loss), ResiduePosition must
    /// round-trip correctly. The limitation is exclusive to custom-loss fragments.
    /// </summary>
    [TestCase(-18.010565, TestName = "NamedLoss_H2O_ResiduePositionPreserved")]
    [TestCase(-17.026549, TestName = "NamedLoss_NH3_ResiduePositionPreserved")]
    [TestCase(-97.976895, TestName = "NamedLoss_H3PO4_ResiduePositionPreserved")]
    [TestCase(0.0, TestName = "NoLoss_ResiduePositionPreserved")]
    public void NamedLoss_ResiduePosition_IsPreservedAfterRoundTrip(double namedLoss)
    {
        const int residuePos = 7;
        string path = TempPath($"named_loss_{System.Math.Abs(namedLoss):F0}");

        MslWriter.Write(path, new[] { EntryWithFragment(namedLoss, residuePos) });

        using MslLibraryData lib = MslReader.Load(path);
        int recovered = lib.Entries[0].MatchedFragmentIons[0].ResiduePosition;

        Assert.That(recovered, Is.EqualTo(residuePos),
            $"ResiduePosition must be preserved for named neutral-loss fragments " +
            $"(NeutralLoss={namedLoss}). The limitation is exclusive to custom losses.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Multiple custom losses in one entry — all lose ResiduePosition
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// When an entry has multiple fragments with different custom losses, every one of
    /// them loses its ResiduePosition individually. The neutral-loss masses themselves
    /// are all recovered correctly from the extended annotation table.
    /// </summary>
    [Test]
    public void CustomLoss_MultipleFragments_AllLoseResiduePosition_MassesStillCorrect()
    {
        const double loss1 = -203.0794; // HexNAc
        const double loss2 = -291.0954; // NeuAc
        string path = TempPath("multi_custom_losses");

        var entry = new MslLibraryEntry
        {
            FullSequence = "PEPTIDE",
            BaseSequence = "PEPTIDE",
            PrecursorMz = 449.74,
            ChargeState = 2,
            MoleculeType = MslFormat.MoleculeType.Peptide,
            DissociationType = DissociationType.HCD,
            MatchedFragmentIons = new List<MslFragmentIon>
            {
                new MslFragmentIon
                {
                    Mz = 300f, Intensity = 1.0f,
                    ProductType = ProductType.y, FragmentNumber = 3,
                    ResiduePosition = 4, Charge = 1, NeutralLoss = loss1
                },
                new MslFragmentIon
                {
                    Mz = 400f, Intensity = 0.8f,
                    ProductType = ProductType.b, FragmentNumber = 5,
                    ResiduePosition = 5, Charge = 1, NeutralLoss = loss2
                }
            }
        };

        MslWriter.Write(path, new[] { entry });

        using MslLibraryData lib = MslReader.Load(path);
        var ions = lib.Entries[0].MatchedFragmentIons;

        // Both must lose ResiduePosition
        Assert.That(ions[0].ResiduePosition, Is.EqualTo(0),
            "First custom-loss fragment must have ResiduePosition=0 after round-trip.");
        Assert.That(ions[1].ResiduePosition, Is.EqualTo(0),
            "Second custom-loss fragment must have ResiduePosition=0 after round-trip.");

        // Both must recover their neutral-loss masses correctly
        var recoveredLosses = new[] { ions[0].NeutralLoss, ions[1].NeutralLoss };
        Assert.That(recoveredLosses, Has.Some.EqualTo(loss1).Within(1e-4),
            $"One of the recovered ions must have neutral loss {loss1}.");
        Assert.That(recoveredLosses, Has.Some.EqualTo(loss2).Within(1e-4),
            $"One of the recovered ions must have neutral loss {loss2}.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Mixed entry: only the custom-loss fragment loses its ResiduePosition
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// In an entry containing both a custom-loss fragment and a named-loss fragment,
    /// only the custom-loss fragment loses its ResiduePosition. The named-loss fragment
    /// is completely unaffected because its ResiduePosition field is never repurposed.
    /// </summary>
    [Test]
    public void MixedEntry_OnlyCustomLossFragment_LosesResiduePosition()
    {
        string path = TempPath("mixed_custom_named");

        var entry = new MslLibraryEntry
        {
            FullSequence = "PEPTIDE",
            BaseSequence = "PEPTIDE",
            PrecursorMz = 449.74,
            ChargeState = 2,
            MoleculeType = MslFormat.MoleculeType.Peptide,
            DissociationType = DissociationType.HCD,
            MatchedFragmentIons = new List<MslFragmentIon>
            {
                // Named loss — ResiduePosition MUST survive
                new MslFragmentIon
                {
                    Mz = 300f, Intensity = 1.0f,
                    ProductType = ProductType.y, FragmentNumber = 3,
                    ResiduePosition = 4, Charge = 1, NeutralLoss = -18.010565
                },
                // Custom loss — ResiduePosition WILL be lost
                new MslFragmentIon
                {
                    Mz = 400f, Intensity = 0.8f,
                    ProductType = ProductType.b, FragmentNumber = 5,
                    ResiduePosition = 5, Charge = 1, NeutralLoss = -203.0794
                }
            }
        };

        MslWriter.Write(path, new[] { entry });

        using MslLibraryData lib = MslReader.Load(path);

        // Identify ions by their recovered neutral-loss mass
        MslFragmentIon? namedIon = null;
        MslFragmentIon? customIon = null;
        foreach (var ion in lib.Entries[0].MatchedFragmentIons)
        {
            if (System.Math.Abs(ion.NeutralLoss - (-18.010565)) < 0.01) namedIon = ion;
            if (System.Math.Abs(ion.NeutralLoss - (-203.0794)) < 0.01) customIon = ion;
        }

        Assert.That(namedIon, Is.Not.Null, "Named-loss ion must be present after round-trip.");
        Assert.That(customIon, Is.Not.Null, "Custom-loss ion must be present after round-trip.");

        Assert.That(namedIon!.ResiduePosition, Is.EqualTo(4),
            "Named-loss fragment ResiduePosition must be preserved " +
            "(not affected by the custom-loss limitation).");
        Assert.That(customIon!.ResiduePosition, Is.EqualTo(0),
            "Custom-loss fragment ResiduePosition must be 0 after round-trip " +
            "(known limitation — field repurposed as ExtAnnotationIdx in binary).");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Writing with ResiduePosition = 0 is a safe no-op
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// If a caller writes a custom-loss fragment with ResiduePosition already set to 0
    /// (i.e. they acknowledge the limitation in advance), the round-trip is clean and the
    /// Debug.Assert does not fire.
    /// </summary>
    [Test]
    public void CustomLoss_ResiduePositionZeroAtWrite_RoundTripClean()
    {
        string path = TempPath("zero_pos_at_write");

        MslWriter.Write(path, new[] { EntryWithFragment(-203.0794, residuePosition: 0) });

        using MslLibraryData lib = MslReader.Load(path);
        var ion = lib.Entries[0].MatchedFragmentIons[0];

        Assert.That(ion.ResiduePosition, Is.EqualTo(0),
            "A custom-loss fragment written with ResiduePosition=0 must read back as 0.");
        Assert.That(ion.NeutralLoss, Is.EqualTo(-203.0794).Within(1e-4),
            "Neutral-loss mass must be correctly recovered when ResiduePosition=0 at write.");
    }
}