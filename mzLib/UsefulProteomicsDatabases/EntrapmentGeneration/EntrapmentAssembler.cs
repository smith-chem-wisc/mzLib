#nullable enable
using MzLibUtil;
using Omics.Digestion;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>What became of one base piece of a target sequence.</summary>
public enum PieceOutcome
{
    /// <summary>Rearranged into an isomeric partner.</summary>
    Permuted,

    /// <summary>
    /// Has no rearrangement, but is shorter than the minimum peptide length, so it is never
    /// identified on its own. Kept exactly as it was: no real peptide enters the entrapment
    /// database, and identity is a rearrangement, so isomerism survives.
    /// </summary>
    KeptVerbatimTooShort,

    /// <summary>
    /// Has no usable rearrangement and *is* long enough to be identified, so it is dropped from the
    /// entrapment sequence. Emitting it unchanged would place a genuine target peptide inside the
    /// entrapment database, where every hit is supposed to be false.
    /// </summary>
    Excised
}

/// <summary>One base piece of the target and what happened to it.</summary>
public sealed class EntrapmentPiece
{
    internal EntrapmentPiece(int index, string targetPiece, string? entrapmentPiece,
        PieceOutcome outcome, EntrapmentFailure failure)
    {
        Index = index;
        TargetPiece = targetPiece;
        EntrapmentPiece_ = entrapmentPiece;
        Outcome = outcome;
        Failure = failure;
    }

    /// <summary>Ordinal of this piece within the target, counted before anything was excised.</summary>
    public int Index { get; }

    public string TargetPiece { get; }

    /// <summary>The piece as it appears in the entrapment sequence, or null when excised.</summary>
    public string? EntrapmentPiece_ { get; }

    public PieceOutcome Outcome { get; }

    /// <summary>Why no partner was found, when <see cref="Outcome"/> is not <see cref="PieceOutcome.Permuted"/>.</summary>
    public EntrapmentFailure Failure { get; }
}

/// <summary>An entrapment sequence, and the record of how it was built.</summary>
public sealed class EntrapmentAssembly
{
    internal EntrapmentAssembly(string targetSequence, string entrapmentSequence,
        int[] targetToEntrapmentPosition, IReadOnlyList<EntrapmentPiece> pieces,
        int missedCleavagePeptidesSpanningAnExcision)
    {
        TargetSequence = targetSequence;
        EntrapmentSequence = entrapmentSequence;
        TargetToEntrapmentPosition = targetToEntrapmentPosition;
        Pieces = pieces;
        MissedCleavagePeptidesSpanningAnExcision = missedCleavagePeptidesSpanningAnExcision;
    }

    public string TargetSequence { get; }

    public string EntrapmentSequence { get; }

    /// <summary>
    /// Position in the entrapment sequence holding each target position's residue, or -1 where the
    /// piece was excised. Modifications ride across on this.
    /// </summary>
    public int[] TargetToEntrapmentPosition { get; }

    public IReadOnlyList<EntrapmentPiece> Pieces { get; }

    /// <summary>
    /// Missed-cleavage peptides of the entrapment sequence whose constituent pieces were not
    /// adjacent in the target, so they have no isomeric counterpart there.
    /// </summary>
    /// <remarks>
    /// Every other missed-cleavage peptide is isomeric for free, because base-piece compositions
    /// add. That is true but not obvious, and these exceptions are counted rather than left
    /// implicit, because anyone computing mass statistics over the entrapment set will assume the
    /// invariant holds everywhere.
    /// </remarks>
    public int MissedCleavagePeptidesSpanningAnExcision { get; }

    public int ExcisedCount => Pieces.Count(p => p.Outcome == PieceOutcome.Excised);

    public int KeptVerbatimCount => Pieces.Count(p => p.Outcome == PieceOutcome.KeptVerbatimTooShort);
}

/// <summary>
/// Builds an entrapment sequence from a target one, by rearranging each of its base pieces in place.
/// </summary>
/// <remarks>
/// <para>The sequence is split at the digestion agent's own cleavage sites
/// (<see cref="DigestionAgent.GetDigestionSiteIndices"/>), each piece is given an isomeric partner
/// by <see cref="EntrapmentPeptideGenerator"/>, and the partners are concatenated back together.
/// Because cleavage sites stay where they were, the entrapment sequence digests into the same
/// pieces, each isomeric with the target's.</para>
/// <para>Rearranging the <b>zero-missed-cleavage</b> pieces is enough to cover missed cleavages
/// too: a missed-cleavage peptide is a run of adjacent pieces, and composition adds, so it is
/// isomeric with its target counterpart without being handled separately. The exception is a run
/// spanning an excision, which has no counterpart -- those are counted.</para>
/// </remarks>
public static class EntrapmentAssembler
{
    /// <summary>
    /// Builds the entrapment sequence for one fold.
    /// </summary>
    /// <param name="targetSequence">The target sequence, typically a protein.</param>
    /// <param name="digestionParams">Supplies the cleavage sites and the minimum peptide length
    /// below which a piece is not identifiable on its own.</param>
    /// <param name="forbiddenSequences">Sequences no partner may equal, normally every target peptide.</param>
    /// <param name="fold">Zero-based fold, in <c>[0, foldCount)</c>.</param>
    /// <param name="foldCount">Partners per target -- the <c>r</c> of an r-fold database.</param>
    /// <param name="seed">Changes every choice reproducibly.</param>
    public static EntrapmentAssembly Assemble(string targetSequence, IDigestionParams digestionParams,
        IReadOnlySet<string> forbiddenSequences, int fold = 0, int foldCount = 1, int seed = 1)
    {
        if (string.IsNullOrEmpty(targetSequence))
        {
            throw new MzLibException("Cannot assemble an entrapment sequence from an empty sequence.");
        }

        List<DigestionMotif> motifs = digestionParams.DigestionAgent.DigestionMotifs;

        // GetDigestionSiteIndices returns 0, every internal cleavage site, and the length -- a
        // complete partition. It comes back from a hash set, so it needs ordering.
        List<int> sites = digestionParams.DigestionAgent.GetDigestionSiteIndices(targetSequence);
        sites.Sort();

        var pieces = new List<EntrapmentPiece>(sites.Count - 1);
        var entrapment = new StringBuilder(targetSequence.Length);
        int[] map = Enumerable.Repeat(-1, targetSequence.Length).ToArray();
        var retainedTargetIndices = new List<int>(sites.Count - 1);

        for (int index = 0; index < sites.Count - 1; index++)
        {
            int start = sites[index];
            int length = sites[index + 1] - start;
            string piece = targetSequence.Substring(start, length);

            EntrapmentPeptide partner = EntrapmentPeptideGenerator.Create(piece, motifs, forbiddenSequences,
                fold, foldCount, seed, TerminalAnchors(index, sites.Count - 1, length));

            if (partner.Succeeded)
            {
                AppendPiece(entrapment, map, start, partner.EntrapmentSequence!, partner.SwappedPositions!);
                pieces.Add(new EntrapmentPiece(index, piece, partner.EntrapmentSequence,
                    PieceOutcome.Permuted, EntrapmentFailure.None));
                retainedTargetIndices.Add(index);
                continue;
            }

            // Too short to be identified on its own: keep it, unchanged, rather than tearing a hole
            // in the protein for a piece nobody could have matched anyway.
            if (piece.Length < digestionParams.MinLength)
            {
                AppendPiece(entrapment, map, start, piece, Identity(piece.Length));
                pieces.Add(new EntrapmentPiece(index, piece, piece,
                    PieceOutcome.KeptVerbatimTooShort, partner.Failure));
                retainedTargetIndices.Add(index);
                continue;
            }

            // Long enough to be identified, and no partner exists: drop it. Its positions stay
            // unmapped, which is what marks them excised.
            pieces.Add(new EntrapmentPiece(index, piece, null, PieceOutcome.Excised, partner.Failure));
        }

        int broken = CountRunsSpanningAGap(retainedTargetIndices, digestionParams.MaxMissedCleavages);

        return new EntrapmentAssembly(targetSequence, entrapment.ToString(), map, pieces, broken);
    }

    /// <summary>
    /// Positions within a piece that must not move because they are the PROTEIN's termini.
    /// </summary>
    /// <remarks>
    /// <para>A modification can be restricted to the N- or C-terminus, and a rearrangement that
    /// moved that residue away would make it invalid for its location -- mzLib then drops it, and
    /// nothing counts the loss. Measured before anchoring: 3,946 N-terminal and 31 C-terminal
    /// modifications lost across the human proteome, 2.4% of all of them.</para>
    /// <para>Both of the first two positions are held, not just the first: a modification annotated
    /// after initiator-methionine cleavage sits on the second residue.</para>
    /// <para>This is unconditional -- the termini are anchored whether or not anything is actually
    /// modified there. Anchoring reactively would make the entrapment sequence a function of
    /// (sequence, seed, modifications) rather than (sequence, seed), so two databases built from the
    /// same proteome with different annotations would disagree on their sequences and the
    /// determinism the pairing rests on would quietly weaken. The cost is two or three positions per
    /// protein, against permutation spaces in the millions.</para>
    /// </remarks>
    private static int[]? TerminalAnchors(int pieceIndex, int pieceCount, int pieceLength)
    {
        bool first = pieceIndex == 0;
        bool last = pieceIndex == pieceCount - 1;
        if (!first && !last)
        {
            return null;
        }

        var anchors = new List<int>(3);
        if (first)
        {
            anchors.Add(0);
            if (pieceLength > 1)
            {
                anchors.Add(1);
            }
        }
        if (last && pieceLength > 0)
        {
            anchors.Add(pieceLength - 1);
        }

        return anchors.ToArray();
    }

    private static void AppendPiece(StringBuilder entrapment, int[] map, int targetStart,
        string entrapmentPiece, int[] swappedWithinPiece)
    {
        int entrapmentStart = entrapment.Length;
        entrapment.Append(entrapmentPiece);

        for (int offset = 0; offset < swappedWithinPiece.Length; offset++)
        {
            map[targetStart + offset] = entrapmentStart + swappedWithinPiece[offset];
        }
    }

    private static int[] Identity(int length) => Enumerable.Range(0, length).ToArray();

    /// <summary>
    /// How many runs of up to <paramref name="maxMissedCleavages"/>+1 retained pieces are not
    /// contiguous in the target -- i.e. how many missed-cleavage peptides straddle an excision.
    /// </summary>
    private static int CountRunsSpanningAGap(List<int> retainedTargetIndices, int maxMissedCleavages)
    {
        int broken = 0;

        for (int i = 0; i < retainedTargetIndices.Count; i++)
        {
            for (int span = 2; span <= maxMissedCleavages + 1 && i + span - 1 < retainedTargetIndices.Count; span++)
            {
                if (retainedTargetIndices[i + span - 1] - retainedTargetIndices[i] != span - 1)
                {
                    broken++;
                }
            }
        }

        return broken;
    }
}
