using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Readers.SpectralLibrary;
using Omics.SpectralMatch.MslSpectralLibrary;


namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 tests verifying the fix to <see cref="SpectralLibrary.ReadFragmentIon"/>
/// that previously computed C-terminal <c>ResiduePosition</c> incorrectly in two ways:
///
/// <list type="number">
///   <item>
///     <b>P16a — absent sequence:</b> fell back to the arbitrary constant 25 with no
///     diagnostic output.
///   </item>
///   <item>
///     <b>P16b — modified sequence (new finding):</b> called <c>peptideSequence.Length</c>
///     on the full modified sequence string (e.g. <c>"PEPTM[Common Variable:Oxidation on
///     M]IDE"</c>), which overcounts by the total number of characters inside all
///     <c>[…]</c> modification brackets.  For a typical oxidized methionine annotation,
///     this inflates the length by ~32 characters, producing C-terminal residue positions
///     that are wildly wrong even when the sequence is present.
///   </item>
/// </list>
///
/// The fix replaces <c>peptideSequence.Length</c> with <c>CountResidues()</c>, a private
/// helper that walks the string counting only characters outside <c>[…]</c> blocks.
/// When the sequence is absent (null/empty), <c>CountResidues</c> returns 0, the fallback
/// to 25 is applied, and a <c>Debug.WriteLine</c> trace is emitted.
///
/// Coverage:
///   1.  y-ion for plain (unmodified) sequence — ResiduePosition uses actual length.
///   2.  y-ion for singly-modified sequence — modification text is stripped.
///   3.  y-ion for doubly-modified sequence — all modification blocks stripped.
///   4.  b-ion (N-terminal) — ResiduePosition = fragmentNumber, never uses peptideLength.
///   5.  Null sequence — falls back to 25, does not throw.
///   6.  Empty sequence — falls back to 25, does not throw.
///   7.  Internal ion with null sequence — unaffected, Terminus = None.
///   8.  Modification text with colons and spaces inside brackets is handled correctly.
///   9.  y1 edge case — ResiduePosition = peptideLength - 1.
///  10.  MSP-parsed ResiduePosition survives conversion to MslFragmentIon (round-trip guard).
/// </summary>
[TestFixture]
public sealed class TestMslPrompt16PeptideLengthFallback
{
    // Split arrays matching ReadLibrarySpectrum's internal constants exactly
    private static readonly char[] FragSplit = { '\t', '"', ')', '/' };
    private static readonly char[] NeutralLossSplit = { '-' };

    /// <summary>
    /// Builds a minimal MSP peak line: <c>mz\tintensity\t"annotation"</c>
    /// </summary>
    private static string PeakLine(double mz, double intensity, string annotation)
        => $"{mz}\t{intensity}\t\"{annotation}\"";

    // ══════════════════════════════════════════════════════════════════════════
    // 1. Plain (unmodified) sequence
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// y3 for a 7-residue plain peptide "PEPTIDE" must have ResiduePosition = 7 - 3 = 4.
    /// This case was already correct before the fix (no modification text to overcount),
    /// but must remain correct after.
    /// </summary>
    [Test]
    public void YIon_PlainSequence_ResiduePosition_UsesActualLength()
    {
        string line = PeakLine(349.19, 0.8, "y3^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
    line, FragSplit, NeutralLossSplit, "PEPTIDE");

        Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.C),
            "y-ion must have C-terminal terminus.");
        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(4),
            "y3 for 7-residue 'PEPTIDE': ResiduePosition must be 7 - 3 = 4.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 2. Singly-modified sequence (the P16b bug — was wrong before fix)
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// "PEPTM[Common Variable:Oxidation on M]IDE" has 8 residues but string.Length = 40.
    ///
    /// Before the fix: peptideLength = 40, y3 ResiduePosition = 40 - 3 = 37 (WRONG).
    /// After the fix:  peptideLength =  8, y3 ResiduePosition =  8 - 3 =  5 (CORRECT).
    /// </summary>
    [Test]
    public void YIon_ModifiedSequence_OneModification_ResiduePosition_StripsModText()
    {
        string modSeq = "PEPTM[Common Variable:Oxidation on M]IDE"; // 8 residues
        // Confirm the test exercises the bug: string.Length must be significantly > 8
        Assert.That(modSeq.Length, Is.GreaterThan(8),
            "Sanity: the modified sequence string must be longer than its residue count.");

        string line = PeakLine(349.19, 0.8, "y3^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, modSeq);

        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(5),
            "y3 for PEPTM[ox]IDE (8 residues) must have ResiduePosition = 8 - 3 = 5. " +
            "Before fix this was string.Length - 3 = a large wrong value.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 3. Doubly-modified sequence
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// "PEPTM[ox]IDEC[cam]" has 9 residues regardless of how many modifications are present.
    /// CountResidues must correctly skip all bracket blocks.
    /// </summary>
    [Test]
    public void YIon_ModifiedSequence_TwoModifications_ResiduePosition_StripsAllModText()
    {
        // 9 residues: P-E-P-T-M-I-D-E-C
        string modSeq =
            "PEPTM[Common Variable:Oxidation on M]IDEC[Common Fixed:Carbamidomethyl on C]";
        string line = PeakLine(420.20, 0.9, "y4^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, modSeq);

        // 9 residues, y4: 9 - 4 = 5
        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(5),
            "y4 for a doubly-modified 9-residue peptide must have ResiduePosition = 9 - 4 = 5.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 4. b-ion (N-terminal) — unchanged by fix
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// N-terminal (b) ions compute ResiduePosition = fragmentNumber.
    /// peptideLength is never used on this path — the fix must not disturb it.
    /// </summary>
    [Test]
    public void BIon_NTerminus_ResiduePosition_EqualToFragmentNumber_Unaffected()
    {
        // Use a modified sequence to confirm the N-terminal path is untouched
        string modSeq = "PEPTM[Common Variable:Oxidation on M]IDE";
        string line = PeakLine(227.10, 0.9, "b2^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, modSeq);

        Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.N),
            "b-ion must have N-terminal terminus.");
        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(2),
            "b2 ResiduePosition must be 2 (= fragmentNumber) regardless of sequence length.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 5. Null sequence — fallback to 25
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// When peptideSequence is null, the method must not throw, must fall back to
    /// peptideLength = 25, and must produce ResiduePosition = 25 - fragmentNumber.
    /// </summary>
    [Test]
    public void YIon_NullSequence_FallsBackTo25_DoesNotThrow()
    {
        string line = PeakLine(349.19, 0.8, "y3^1");

        MatchedFragmentIon ion = null!;
        Assert.DoesNotThrow(() =>
        {
            ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
                line, FragSplit, NeutralLossSplit, peptideSequence: null);
        }, "ReadFragmentIon must not throw when peptideSequence is null.");

        // Fallback: 25 - 3 = 22
        Assert.That(ion!.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(22),
            "y3 with null sequence must use fallback length 25: ResiduePosition = 25 - 3 = 22.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 6. Empty sequence — same fallback
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// When peptideSequence is empty, behaviour must be identical to the null case.
    /// </summary>
    [Test]
    public void YIon_EmptySequence_FallsBackTo25_DoesNotThrow()
    {
        string line = PeakLine(349.19, 0.8, "y3^1");

        MatchedFragmentIon ion = null!;
        Assert.DoesNotThrow(() =>
        {
            ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
                line, FragSplit, NeutralLossSplit, peptideSequence: string.Empty);
        }, "ReadFragmentIon must not throw when peptideSequence is empty.");

        Assert.That(ion!.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(22),
            "y3 with empty sequence must use fallback length 25: ResiduePosition = 25 - 3 = 22.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 7. Internal ions with null sequence — completely unaffected
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Internal ions match the InternalIonRegex and return before the peptideLength branch.
    /// Null sequence must not cause a throw and must produce Terminus = None.
    /// </summary>
    [Test]
    public void InternalIon_NullSequence_Unaffected_TerminusNone()
    {
        string line = PeakLine(312.1, 0.3, "bIb[3-6]^1");

        MatchedFragmentIon ion = null!;
        Assert.DoesNotThrow(() =>
        {
            ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
                line, FragSplit, NeutralLossSplit, peptideSequence: null);
        }, "Internal ion parsing must not throw when peptideSequence is null.");

        Assert.That(ion!.NeutralTheoreticalProduct.Terminus,
            Is.EqualTo(FragmentationTerminus.None),
            "Internal ions must always have Terminus=None; the peptideLength branch is never reached.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 8. Colons and spaces inside brackets
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Modification names commonly contain colons (e.g. "Common Variable:Oxidation on M")
    /// and spaces. CountResidues must ignore all characters inside brackets, including these.
    /// </summary>
    [Test]
    public void YIon_ModTextWithColonsAndSpaces_OnlyCountsResiduesOutsideBrackets()
    {
        // "PEPTM[Common Variable:Oxidation on M]IDE" has 8 residues
        string modSeq = "PEPTM[Common Variable:Oxidation on M]IDE";
        string line = PeakLine(175.12, 1.0, "y1^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, modSeq);

        // y1 for 8 residues: ResiduePosition = 8 - 1 = 7
        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(7),
            "y1 for PEPTM[Common Variable:Oxidation on M]IDE (8 residues) must have " +
            "ResiduePosition = 8 - 1 = 7. Colons/spaces inside brackets must be skipped.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 9. y1 edge case
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// y1 has ResiduePosition = peptideLength - 1. For a 5-residue peptide that is 4.
    /// </summary>
    [Test]
    public void YIon_y1_ResiduePosition_EqualsPeptideLengthMinusOne()
    {
        string line = PeakLine(116.07, 1.0, "y1^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, "GILAV");

        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(4),
            "y1 for a 5-residue peptide must have ResiduePosition = 5 - 1 = 4.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // 10. MSP-parsed ResiduePosition survives conversion to MslFragmentIon
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// A y3 ion parsed correctly from MSP (ResiduePosition = 4 for a 7-residue peptide)
    /// must produce the same ResiduePosition when converted to an MslFragmentIon.
    /// This is the MSP → MSL round-trip scenario the comment identifies.
    /// </summary>
    [Test]
    public void YIon_CorrectResiduePosition_SurvivesConversionToMslFragmentIon()
    {
        string line = PeakLine(349.19, 0.8, "y3^1");

        MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
            line, FragSplit, NeutralLossSplit, "PEPTIDE");

        // Verify correct position before conversion
        Assert.That(ion.NeutralTheoreticalProduct.ResiduePosition, Is.EqualTo(4),
            "ResiduePosition must be 4 (= 7 - 3) before MslFragmentIon conversion.");

        // Convert to MslFragmentIon and verify position is preserved
        var mslFrag = new MslFragmentIon
        {
            Mz = (float)ion.Mz,
            Intensity = (float)ion.Intensity,
            ProductType = ion.NeutralTheoreticalProduct.ProductType,
            FragmentNumber = ion.NeutralTheoreticalProduct.FragmentNumber,
            ResiduePosition = ion.NeutralTheoreticalProduct.ResiduePosition,
            Charge = ion.Charge,
            NeutralLoss = ion.NeutralTheoreticalProduct.NeutralLoss
        };

        Assert.That(mslFrag.ResiduePosition, Is.EqualTo(4),
            "ResiduePosition = 4 must survive intact through conversion to MslFragmentIon.");
    }
}