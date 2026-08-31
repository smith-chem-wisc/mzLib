#nullable enable
using MzLibUtil;
using Omics.Digestion;
using System;
using System.Collections.Generic;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>
/// Why an entrapment peptide could not be produced. These are kept apart because the caller's
/// remedy differs for each, and merging them into a single "failed" hides that.
/// </summary>
public enum EntrapmentFailure
{
    /// <summary>A partner was produced.</summary>
    None,

    /// <summary>
    /// The sequence has exactly one arrangement once its cleavage sites are held in place -- a
    /// homopolymeric tract such as "SSSSSSR". Arithmetic, not competition: no seed, no extra
    /// searching and no smaller fold count can produce a partner for this peptide.
    /// </summary>
    NoPermutationExists,

    /// <summary>
    /// Arrangements exist, but every one this fold may draw on is already spoken for -- it is a
    /// target peptide, or it has already been issued. Depends on the database, so a different
    /// target set may succeed where this one did not.
    /// </summary>
    AllPermutationsTaken,

    /// <summary>
    /// The space is real but too small to give every fold its own partner. The remedy is a smaller
    /// fold count, which distinguishes this from <see cref="AllPermutationsTaken"/>.
    /// </summary>
    SpaceTooSmallForFoldCount
}

/// <summary>
/// One target peptide's entrapment partner, or the reason it has none.
/// </summary>
public sealed class EntrapmentPeptide
{
    internal EntrapmentPeptide(string targetSequence, string? entrapmentSequence, int[]? swappedPositions,
        int fold, BigInteger permutationSpaceSize, int probesUsed, EntrapmentFailure failure)
    {
        TargetSequence = targetSequence;
        EntrapmentSequence = entrapmentSequence;
        SwappedPositions = swappedPositions;
        Fold = fold;
        PermutationSpaceSize = permutationSpaceSize;
        ProbesUsed = probesUsed;
        Failure = failure;
    }

    public string TargetSequence { get; }

    /// <summary>The partner, isomeric with the target. Null when <see cref="Succeeded"/> is false.</summary>
    public string? EntrapmentSequence { get; }

    /// <summary>
    /// Maps each position in the target to its position in the partner, the same contract as
    /// <see cref="DecoySequenceValidator.ScrambleSequence"/>, so modifications ride across unchanged.
    /// </summary>
    public int[]? SwappedPositions { get; }

    /// <summary>Zero-based fold this partner belongs to.</summary>
    public int Fold { get; }

    /// <summary>How many distinct arrangements exist. Reported even on failure, because the value
    /// is what distinguishes an impossible peptide from a merely crowded one.</summary>
    public BigInteger PermutationSpaceSize { get; }

    /// <summary>Candidates examined, including the one returned. One means the first choice was free.</summary>
    public int ProbesUsed { get; }

    public EntrapmentFailure Failure { get; }

    public bool Succeeded => Failure == EntrapmentFailure.None;
}

/// <summary>
/// Produces the entrapment partner of a target peptide: the same residues in a different order,
/// with every cleavage site left where it was.
/// </summary>
/// <remarks>
/// <para>The partner is <b>isomeric</b> with its target -- same composition, so same mass and the
/// same count of any residue of interest. For localization work that last property is the point:
/// the number of candidate sites in a peptide sets the difficulty of placing a modification on it,
/// so a partner carrying a different count would be compared at a difficulty its target never had.
/// It also makes target and partner the hardest possible pair to tell apart, which is the
/// discrimination an entrapment experiment is trying to measure.</para>
/// <para>Selection is by <b>index, not by chance</b>. The starting point is derived from the
/// sequence and the seed, and search walks forward from there, so the answer is a pure function of
/// its inputs: unchanged by processing order, by the protein the peptide came from, by threading,
/// and by framework version. Nothing here consumes random state, and folds never consult one
/// another, so a database can be regenerated in pieces without invalidating what already exists.
/// </para>
/// <para>Because the search enumerates rather than samples, "no partner exists" is a
/// <b>proven</b> statement rather than a search that gave up, and the three ways it can fail are
/// reported apart (<see cref="EntrapmentFailure"/>).</para>
/// </remarks>
public static class EntrapmentPeptideGenerator
{
    /// <summary>
    /// The entrapment partner of <paramref name="targetSequence"/> for one fold.
    /// </summary>
    /// <param name="targetSequence">The target peptide.</param>
    /// <param name="motifs">Cleavage motifs whose residues must stay in place, so that the partner
    /// still digests the same way. Typically <c>DigestionAgent.DigestionMotifs</c>.</param>
    /// <param name="forbiddenSequences">Sequences the partner may not equal -- normally every target
    /// peptide in the database.</param>
    /// <param name="fold">Zero-based fold, in <c>[0, foldCount)</c>.</param>
    /// <param name="foldCount">How many partners each target is to receive (the <c>r</c> of an
    /// r-fold entrapment database).</param>
    /// <param name="seed">Changes every choice reproducibly.</param>
    /// <param name="alsoHeldInPlace">Extra zero-based positions within this peptide to hold still
    /// on top of the cleavage sites. The assembler uses this to anchor a protein's termini, so a
    /// modification restricted to one of them survives the rearrangement.</param>
    public static EntrapmentPeptide Create(string targetSequence, List<DigestionMotif> motifs,
        IReadOnlySet<string> forbiddenSequences, int fold = 0, int foldCount = 1, int seed = 1,
        IReadOnlyCollection<int>? alsoHeldInPlace = null)
    {
        if (string.IsNullOrEmpty(targetSequence))
        {
            throw new MzLibException("Cannot build an entrapment peptide from an empty sequence.");
        }
        if (foldCount < 1)
        {
            throw new MzLibException($"Fold count must be at least 1, but was {foldCount}.");
        }
        if (fold < 0 || fold >= foldCount)
        {
            throw new MzLibException($"Fold {fold} is outside the {foldCount} requested folds.");
        }

        BigInteger size = DecoySequenceValidator.PermutationSpaceSize(targetSequence, motifs, alsoHeldInPlace);

        // One arrangement means the identity and nothing else. No fold count and no seed can help.
        if (size <= BigInteger.One)
        {
            return Failed(targetSequence, fold, size, EntrapmentFailure.NoPermutationExists);
        }

        // Give each fold its own contiguous stretch of the space. Disjoint stretches make the folds
        // distinct by construction and independent of one another -- neither has to know what the
        // others chose, so they can be produced in any order, in parallel, or years apart.
        BigInteger stretch = size / foldCount;
        if (stretch.IsZero)
        {
            return Failed(targetSequence, fold, size, EntrapmentFailure.SpaceTooSmallForFoldCount);
        }

        BigInteger start = fold * stretch;
        BigInteger offset = DeriveOffset(targetSequence, seed, stretch);

        for (BigInteger step = BigInteger.Zero; step < stretch; step++)
        {
            BigInteger index = start + (offset + step) % stretch;
            string candidate = DecoySequenceValidator.UnrankPermutation(targetSequence, motifs, index,
                out int[] swapped, alsoHeldInPlace);

            if (candidate != targetSequence && !forbiddenSequences.Contains(candidate))
            {
                return new EntrapmentPeptide(targetSequence, candidate, swapped, fold, size,
                    (int)(step + BigInteger.One), EntrapmentFailure.None);
            }
        }

        // The stretch was walked end to end, so this is a proof rather than an abandoned search.
        return Failed(targetSequence, fold, size, EntrapmentFailure.AllPermutationsTaken);
    }

    private static EntrapmentPeptide Failed(string targetSequence, int fold, BigInteger size,
        EntrapmentFailure failure) =>
        new(targetSequence, null, null, fold, size, 0, failure);

    /// <summary>
    /// Where in a fold's stretch to start looking, derived from the sequence and the seed.
    /// </summary>
    /// <remarks>
    /// SHA-256 rather than <see cref="string.GetHashCode()"/> or a home-made hash: its output is
    /// fixed by specification, so a database regenerated on another machine, another runtime or in
    /// another decade is byte-identical. String hash codes are explicitly not stable across runs.
    /// </remarks>
    private static BigInteger DeriveOffset(string sequence, int seed, BigInteger stretch)
    {
        byte[] digest = SHA256.HashData(Encoding.UTF8.GetBytes($"{seed}:{sequence}"));

        // Leading zero byte keeps BigInteger's two's-complement reading non-negative.
        byte[] unsigned = new byte[digest.Length + 1];
        Array.Copy(digest, unsigned, digest.Length);

        return new BigInteger(unsigned) % stretch;
    }
}
