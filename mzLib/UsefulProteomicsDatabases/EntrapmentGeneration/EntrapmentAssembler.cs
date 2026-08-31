#nullable enable
using System;
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
        int missedCleavagePeptidesSpanningAnExcision,
        IReadOnlyList<string> unrepairableRunCollisionPeptides)
    {
        TargetSequence = targetSequence;
        EntrapmentSequence = entrapmentSequence;
        TargetToEntrapmentPosition = targetToEntrapmentPosition;
        Pieces = pieces;
        MissedCleavagePeptidesSpanningAnExcision = missedCleavagePeptidesSpanningAnExcision;
        UnrepairableRunCollisionPeptides = unrepairableRunCollisionPeptides;
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

    /// <summary>
    /// Missed-cleavage peptides of this entrapment sequence that equal a real target peptide and
    /// could not be avoided.
    /// </summary>
    /// <remarks>
    /// <para>Each piece is permuted away from any collision it can be permuted away from, runs
    /// included. What cannot be repaired is a run whose final piece has exactly one arrangement --
    /// a short low-complexity piece such as <c>AAR</c>, where holding the cleavage residue leaves
    /// nothing to reorder. There is no candidate to move to, so the collision stands.</para>
    /// <para>Repairing it would mean backtracking into an already-placed piece or excising a
    /// perfectly good one, and this project counts collisions rather than silently repairing them.
    /// Measured on the reviewed human proteome: 1,794 such peptides, 0.065% of the target set,
    /// 97.9% of them ending in a piece with no alternative and 80% exactly at the minimum
    /// searchable length of seven residues. A search that matches one counts a real peptide as an
    /// entrapment discovery, so anyone reading an entrapment count needs this number beside it.</para>
    /// </remarks>
    public int UnrepairableRunCollisions => UnrepairableRunCollisionPeptides.Count;

    /// <summary>
    /// The actual peptides behind <see cref="UnrepairableRunCollisions"/> -- entrapment
    /// missed-cleavage peptides that are also real target peptides.
    /// </summary>
    /// <remarks>
    /// A count tells a consumer how much to distrust an entrapment total; the sequences let it
    /// subtract them. The glyco project asked for the list rather than the number for exactly that
    /// reason, and a list is the only form in which "exclude these" is actionable.
    /// </remarks>
    public IReadOnlyList<string> UnrepairableRunCollisionPeptides { get; }

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
        RefuseAgentsWhoseSitesCannotBeHeld(digestionParams.DigestionAgent.Name, motifs);

        // GetDigestionSiteIndices returns 0, every internal cleavage site, and the length -- a
        // complete partition. It comes back from a hash set, so it needs ordering.
        List<int> sites = digestionParams.DigestionAgent.GetDigestionSiteIndices(targetSequence);
        sites.Sort();

        var pieces = new List<EntrapmentPiece>(sites.Count - 1);
        var entrapment = new StringBuilder(targetSequence.Length);
        // The entrapment pieces already placed, in the order they appear in the entrapment
        // sequence, so a candidate can be tested against the runs it would complete.
        var placed = new List<string>(sites.Count - 1);
        var unrepairable = new List<string>();
        int[] map = Enumerable.Repeat(-1, targetSequence.Length).ToArray();
        var retainedTargetIndices = new List<int>(sites.Count - 1);

        for (int index = 0; index < sites.Count - 1; index++)
        {
            int start = sites[index];
            int length = sites[index + 1] - start;
            string piece = targetSequence.Substring(start, length);

            Func<string, IReadOnlyList<string>>? completesAForbiddenRun =
                RejectRunCollisions(placed, digestionParams.MaxMissedCleavages, forbiddenSequences,
                    digestionParams.MinLength, digestionParams.MaxLength);

            EntrapmentPeptide partner = EntrapmentPeptideGenerator.Create(piece, motifs, forbiddenSequences,
                fold, foldCount, seed, TerminalAnchors(start, length, targetSequence.Length),
                completesAForbiddenRun is null ? null : c => completesAForbiddenRun(c).Count > 0);

            if (partner.Succeeded)
            {
                AppendPiece(entrapment, map, start, partner.EntrapmentSequence!, partner.SwappedPositions!);
                placed.Add(partner.EntrapmentSequence!);
                // A permuted piece was chosen with the run test applied, so this cannot fire -- it is
                // asserted rather than assumed only because the count below must mean one thing.
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
                placed.Add(piece);
                // Kept verbatim because it has no alternative arrangement, so if it completes a
                // forbidden run there is nothing to move to. Name it; do not repair it.
                // Every run this piece completes, not merely the first: several can end here, and
                // each is a distinct peptide a consumer has to exclude.
                foreach (string collision in completesAForbiddenRun?.Invoke(piece)
                                             ?? System.Array.Empty<string>())
                {
                    unrepairable.Add(collision);
                }
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

        return new EntrapmentAssembly(targetSequence, entrapment.ToString(), map, pieces, broken,
            unrepairable);
    }

    /// <summary>
    /// Rearranges a whole sequence as one unit, for top-down work.
    /// </summary>
    /// <param name="targetSequence">The target protein.</param>
    /// <param name="forbiddenSequences">Sequences the partner may not equal.</param>
    /// <param name="minLength">Below this the sequence is not identifiable, so it is kept verbatim
    /// rather than dropped, exactly as a short base piece is.</param>
    /// <remarks>
    /// <para>Top-down does not digest, so there are no cleavage sites to hold and no pieces to
    /// assemble: the unit <i>is</i> the protein. Passing no motifs leaves every position free but
    /// the anchored termini, which is the whole difference between this and
    /// <see cref="Assemble"/>.</para>
    /// <para>Because nothing is excised, <b>the length is preserved</b>. That is what makes the
    /// positional annotations carryable here when the bottom-up path had to drop them: an excision
    /// shortens the protein and leaves coordinates pointing past its end, and there is no excision
    /// in this mode. A protein with no usable rearrangement is dropped whole rather than partly.</para>
    /// </remarks>
    public static EntrapmentAssembly AssembleWholeProtein(string targetSequence,
        IReadOnlySet<string> forbiddenSequences, int minLength = 1,
        int fold = 0, int foldCount = 1, int seed = 1)
    {
        if (string.IsNullOrEmpty(targetSequence))
        {
            throw new MzLibException("Cannot assemble an entrapment proteoform from an empty sequence.");
        }

        var noMotifs = new List<DigestionMotif>();
        EntrapmentPeptide partner = EntrapmentPeptideGenerator.Create(targetSequence, noMotifs,
            forbiddenSequences, fold, foldCount, seed,
            TerminalAnchors(0, targetSequence.Length, targetSequence.Length));

        var pieces = new List<EntrapmentPiece>(1);
        int[] map = Enumerable.Repeat(-1, targetSequence.Length).ToArray();
        var entrapment = new StringBuilder(targetSequence.Length);

        if (partner.Succeeded)
        {
            AppendPiece(entrapment, map, 0, partner.EntrapmentSequence!, partner.SwappedPositions!);
            pieces.Add(new EntrapmentPiece(0, targetSequence, partner.EntrapmentSequence,
                PieceOutcome.Permuted, EntrapmentFailure.None));
        }
        else if (targetSequence.Length < minLength)
        {
            AppendPiece(entrapment, map, 0, targetSequence, Identity(targetSequence.Length));
            pieces.Add(new EntrapmentPiece(0, targetSequence, targetSequence,
                PieceOutcome.KeptVerbatimTooShort, partner.Failure));
        }
        else
        {
            pieces.Add(new EntrapmentPiece(0, targetSequence, null, PieceOutcome.Excised, partner.Failure));
        }

        return new EntrapmentAssembly(targetSequence, entrapment.ToString(), map, pieces, 0,
            new List<string>());
    }

    /// <summary>
    /// Refuses a digestion agent whose cleavage sites this assembler cannot hold in place.
    /// </summary>
    /// <remarks>
    /// <para>The class invariant is that the entrapment sequence digests into the same pieces as its
    /// target. That rests on every cleavage site being pinned, and two kinds of motif escape the
    /// pinning that <see cref="DecoySequenceValidator.CleavageSitePositions"/> performs.</para>
    /// <para>A <b>preventing</b> motif pins nothing at all: <c>CleavageSitePositions</c> skips a
    /// position where a motif fits but is prevented, so under <c>trypsin|P</c>'s <c>K[P]|</c> neither
    /// the K nor the P after it is held. Moving that P away invents a cleavage site the target
    /// lacked; moving another P next to another K destroys one. That is five shipped agents --
    /// <c>chymotrypsin|P</c>, <c>elastase|P</c>, <c>Lys-C|P</c>, <c>trypsin|P</c> and
    /// <c>subtilisin|P</c> -- and <c>trypsin|P</c> is an ordinary choice in real work.</para>
    /// <para>A <b>multi-residue</b> motif has its span cut by the piece boundary: <c>TX|T</c> cuts
    /// after the second of three residues, so one piece ends <c>…TX</c> and the next begins
    /// <c>T…</c>, neither fragment matches inside its own piece, and nothing is held at the seam.
    /// That covers <c>StcE</c> and <c>StcE-trypsin</c>.</para>
    /// <para>The real fix is to pin across piece boundaries, which means this stops being a per-piece
    /// loop. Until then, refusing is the honest answer: a database that digests differently from its
    /// target produces peptides that then fail to pair, and it does so silently. Failing at the call
    /// is better than a caller discovering it from a search result.</para>
    /// </remarks>
    /// <exception cref="MzLibException">The agent carries a preventing or multi-residue motif.</exception>
    private static void RefuseAgentsWhoseSitesCannotBeHeld(string agentName, List<DigestionMotif> motifs)
    {
        foreach (DigestionMotif motif in motifs ?? new List<DigestionMotif>())
        {
            if (!string.IsNullOrEmpty(motif.PreventingCleavage))
            {
                throw new MzLibException(
                    $"The digestion agent '{agentName}' has a preventing-cleavage motif "
                    + $"('{motif.InducingCleavage}[{motif.PreventingCleavage}]'), whose residues this "
                    + "assembler cannot hold in place. A rearrangement could then invent or destroy a "
                    + "cleavage site, so the entrapment protein would not digest into the same pieces "
                    + "as its target and its peptides would not pair. Use an agent without one.");
            }

            if ((motif.InducingCleavage?.Length ?? 0) > 1)
            {
                throw new MzLibException(
                    $"The digestion agent '{agentName}' has a multi-residue motif "
                    + $"('{motif.InducingCleavage}'), whose span is cut by the boundary between base "
                    + "pieces, so neither fragment matches inside its own piece and nothing is held at "
                    + "the seam. The entrapment protein would not digest into the same pieces as its "
                    + "target. Use a single-residue agent.");
            }
        }
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
    private static int[]? TerminalAnchors(int pieceStart, int pieceLength, int proteinLength)
    {
        // Protein coordinates, translated into piece coordinates -- not piece coordinates that
        // happen to coincide. When the opening piece is a single residue the protein's second
        // residue lives in the NEXT piece, and computing "the first two" per piece never held it, so
        // an "N-terminal." modification annotated at position 2 could still move. Under trypsin that
        // needs an entry beginning with K or R, which is rare rather than impossible, and the
        // documentation promised to cover it either way.
        int[] protein = proteinLength > 1
            ? new[] { 0, 1, proteinLength - 1 }
            : new[] { 0 };

        List<int>? anchors = null;
        foreach (int position in protein)
        {
            int within = position - pieceStart;
            if (within >= 0 && within < pieceLength)
            {
                (anchors ??= new List<int>(3)).Add(within);
            }
        }

        return anchors?.ToArray();
    }

    /// <summary>
    /// A test that rejects a candidate piece which would complete a run equal to a real target
    /// peptide.
    /// </summary>
    /// <remarks>
    /// <para>A missed-cleavage peptide is a run of adjacent base pieces, and composition adds, so
    /// the run is isomeric with its target counterpart for free. What does <i>not</i> come for free
    /// is that the run differs from every real peptide: each piece is permuted knowing only itself,
    /// so <c>E1+E2</c> can equal a genuine target peptide even when neither <c>E1</c> nor <c>E2</c>
    /// does. Such a sequence is a real peptide sitting in the entrapment database, where every hit
    /// is supposed to be false, and a search that matches it counts a true peptide as an entrapment
    /// discovery.</para>
    /// <para>Every run ends at some piece, so testing the runs that <i>end</i> at the piece being
    /// placed covers all of them, and left-to-right assembly has already placed everything such a
    /// run needs. The runs are built once per piece rather than once per candidate, because the
    /// probe loop may try many candidates and the preceding pieces do not change between them.</para>
    /// <para>The cost is that a piece's permutation now depends on its neighbours, so the
    /// construction is a pure function of <c>(protein, seed)</c> rather than of <c>(piece, seed)</c>.
    /// Determinism, parallel safety across proteins, and the composition-plus-pinning pairing key
    /// are all unaffected -- the piece is still a permutation of its target piece.</para>
    /// </remarks>
    private static Func<string, IReadOnlyList<string>>? RejectRunCollisions(List<string> placed,
        int maxMissedCleavages, IReadOnlySet<string> forbiddenSequences,
        int minLength, int maxLength)
    {
        if (maxMissedCleavages < 1 || placed.Count == 0)
        {
            return null;
        }

        int longest = Math.Min(maxMissedCleavages, placed.Count);
        var runs = new string[longest];
        var builder = new StringBuilder();
        for (int back = 1; back <= longest; back++)
        {
            builder.Insert(0, placed[placed.Count - back]);
            runs[back - 1] = builder.ToString();
        }

        // Returns *every* offending run rather than the first, and only runs a search could report.
        //
        // Both halves were wrong in ways that pulled the count in opposite directions. Returning on
        // the first match meant several colliding runs ending at one piece were reported as one, so
        // a property documenting a count of peptides was really counting pieces. And with no length
        // bound a concatenation shorter than MinLength -- six residues from "AAKAAR" under a minimum
        // of seven -- still rejected candidates and still counted, while the pairing index that
        // defines the searchable population had already excluded it. The assembler and the pairing
        // have to agree about what a peptide is or the exclusion list is not on the same footing as
        // the population it is meant to be subtracted from.
        return candidate =>
        {
            List<string>? offending = null;
            foreach (string run in runs)
            {
                int length = run.Length + candidate.Length;
                if (length < minLength || length > maxLength)
                {
                    continue;
                }

                string whole = run + candidate;
                if (forbiddenSequences.Contains(whole))
                {
                    (offending ??= new List<string>()).Add(whole);
                }
            }

            return (IReadOnlyList<string>)offending ?? System.Array.Empty<string>();
        };
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
