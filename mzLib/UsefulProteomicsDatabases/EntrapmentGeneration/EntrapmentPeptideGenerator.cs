#nullable enable
using MzLibUtil;
using Omics.Digestion;
using System;
using System.Globalization;
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
    /// Every arrangement in this fold's share was refused for something outside the piece itself:
    /// it would have made a <b>missed-cleavage</b> peptide equal to a real target peptide, or -- for
    /// the piece that opens a protein -- its <b>initiator-methionine-stripped</b> form would have
    /// been one.
    /// </summary>
    /// <remarks>
    /// Kept apart from <see cref="AllPermutationsTaken"/> because the remedies differ, which is the
    /// whole reason these causes are separate. A peptide whose own arrangements are spoken for wants
    /// a different target database; a peptide defeated by the runs it would complete wants a
    /// different seed, or different neighbours, and is a property of where it sits rather than of
    /// what it is.
    /// </remarks>
    RunCollisionsExhaustedTheSpace,

    /// <summary>
    /// The sequence has exactly one arrangement once its cleavage sites are held in place -- a
    /// homopolymeric tract such as "SSSSSSR". Arithmetic, not competition: no seed, no extra
    /// searching and no smaller fold count can produce a partner for this peptide.
    /// </summary>
    NoPermutationExists,

    /// <summary>
    /// Arrangements exist, but every one this fold may draw on is a forbidden sequence -- normally a
    /// target peptide. Nothing tracks partners already issued: folds draw on disjoint shares, so they
    /// cannot take one another's. Depends on the database, so a different target set may succeed
    /// where this one did not.
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

    /// <summary>
    /// Candidates examined, including the one returned. One means the first choice was free; on a
    /// failure it is the whole share, which is what makes the failure a proof. Saturates at
    /// <see cref="int.MaxValue"/> rather than overflowing, being a diagnostic.
    /// </summary>
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
/// <para>Selection is by <b>index, not by chance</b>. Where a fold starts in its share, and the
/// keyed order it walks the share in, are both derived from the sequence and the seed, so the answer
/// is a pure function of the arguments: unchanged by processing order, by threading, and by
/// framework version. The held positions and the context test are arguments too, and they are where
/// the surrounding protein enters -- a piece that opens a protein is held differently from the same
/// piece inside one, so the two can be given different partners. Nothing here consumes random
/// state, and folds never consult one another, so a database can be regenerated in pieces without
/// invalidating what already exists.</para>
/// <para>Because the search enumerates rather than samples, "no partner exists" is a
/// <b>proven</b> statement rather than a search that gave up, and the four ways it can fail are
/// reported apart (<see cref="EntrapmentFailure"/>).</para>
/// </remarks>
public static class EntrapmentPeptideGenerator
{
    /// <summary>Stands in for a null exclusion set, so "nothing is forbidden" costs no allocation.</summary>
    private static readonly IReadOnlySet<string> NoForbiddenSequences = new HashSet<string>();

    /// <summary>
    /// The entrapment partner of <paramref name="targetSequence"/> for one fold.
    /// </summary>
    /// <param name="targetSequence">The target peptide.</param>
    /// <param name="motifs">Cleavage motifs whose residues must stay in place, so that the partner
    /// still digests the same way. Typically <c>DigestionAgent.DigestionMotifs</c>.</param>
    /// <param name="forbiddenSequences">Sequences the partner may not equal -- normally every target
    /// peptide in the database. Null is read as "nothing is forbidden", the same as an empty set.</param>
    /// <param name="fold">Zero-based fold, in <c>[0, foldCount)</c>.</param>
    /// <param name="foldCount">How many partners each target is to receive (the <c>r</c> of an
    /// r-fold entrapment database).</param>
    /// <param name="seed">Changes every choice reproducibly.</param>
    /// <param name="alsoHeldInPlace">Extra zero-based positions within this peptide to hold still
    /// on top of the cleavage sites. The assembler uses this to anchor a protein's termini, so a
    /// modification restricted to one of them survives the rearrangement.</param>
    /// <param name="rejectInContext">An extra test a candidate must pass, on top of not being a
    /// forbidden sequence itself. The assembler uses this to reject a candidate that would make a
    /// *missed-cleavage* peptide equal to a real target peptide: this generator sees one base piece
    /// at a time, and a missed-cleavage peptide is a run of adjacent pieces, so a concatenation can
    /// collide even when none of its parts does. Measured before this existed: 2,094 peptides
    /// (0.076%) appeared in both the target and entrapment sets, 2,090 of them with a missed
    /// cleavage.</param>
    /// <exception cref="MzLibException">The motifs are missing, or include one whose sites cannot be
    /// held in place (empty, preventing, multi-residue or wildcard), so a partner could digest
    /// differently from its target; or an argument is out of range.</exception>
    public static EntrapmentPeptide Create(string targetSequence, List<DigestionMotif> motifs,
        IReadOnlySet<string> forbiddenSequences, int fold = 0, int foldCount = 1, int seed = 1,
        IReadOnlyCollection<int>? alsoHeldInPlace = null, Func<string, bool>? rejectInContext = null)
    {
        // Vetted here as well as by the assembler, because this is public and the class promises
        // every cleavage site stays where it was. Under trypsin|P the vetting lived only on the
        // protein path, so a caller of this method got AKLPPR -> AKPPLR, a partner missing one of
        // its target's two cleavage sites.
        if (motifs is null)
        {
            throw new MzLibException(
                "Cannot hold cleavage sites in place without the cleavage motifs; pass the digestion "
                + "agent's DigestionMotifs.");
        }
        EntrapmentAssembler.RefuseAgentsWhoseSitesCannotBeHeld(null, motifs);

        return CreateFromVettedMotifs(targetSequence, motifs, forbiddenSequences, fold, foldCount, seed,
            alsoHeldInPlace, rejectInContext);
    }

    /// <summary>
    /// <see cref="Create"/> without vetting the motifs, for a caller that has already refused agents
    /// whose sites cannot be held: the assembler, which vets once per protein rather than once per
    /// piece and fold.
    /// </summary>
    internal static EntrapmentPeptide CreateFromVettedMotifs(string targetSequence, List<DigestionMotif> motifs,
        IReadOnlySet<string> forbiddenSequences, int fold, int foldCount, int seed,
        IReadOnlyCollection<int>? alsoHeldInPlace, Func<string, bool>? rejectInContext)
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

        // No exclusions rather than an exception. The parameter is optional in spirit -- a caller
        // generating against no target database has nothing to forbid -- and reaching `.Contains`
        // on a null set threw a NullReferenceException out of a public API, which says nothing
        // about which argument was wrong.
        forbiddenSequences ??= NoForbiddenSequences;

        BigInteger identity = DecoySequenceValidator.RankPermutation(targetSequence, motifs, alsoHeldInPlace,
            out BigInteger size);

        // One arrangement means the identity and nothing else. No fold count and no seed can help.
        if (size <= BigInteger.One)
        {
            return Failed(targetSequence, fold, size, 0, EntrapmentFailure.NoPermutationExists);
        }

        // The identity is never a usable partner, so the space a fold count has to be shared out
        // over is `size - 1`, not `size`. Testing `size / foldCount == 0` missed the boundary: at
        // size 3 and foldCount 3 every fold gets a stretch of one, and whichever stretch holds the
        // identity has that single candidate refused and reported as AllPermutationsTaken -- which
        // sends a caller after a different target database when the answer is a smaller fold count.
        if (size - BigInteger.One < foldCount)
        {
            return Failed(targetSequence, fold, size, 0, EntrapmentFailure.SpaceTooSmallForFoldCount);
        }

        // Share out only the arrangements that differ from the target, and share them by residue
        // class modulo the fold count: fold f owns non-identity ranks f, f + foldCount, f + 2·foldCount,
        // ... Disjoint shares make the folds distinct by construction and independent of one
        // another -- neither has to know what the others chose, so they can be produced in any
        // order, in parallel, or years apart.
        //
        // Interleaved rather than contiguous (Alexander-Sol, #1271): `size / foldCount` left the
        // remainder unused while the identity still sat inside one fold's block, so at a block of
        // one that fold was excised as SpaceTooSmallForFoldCount although unused arrangements
        // remained. Skipping the identity's rank and spreading the remainder over the low folds
        // means every fold holds at least one candidate whenever the guard above passes.
        //
        // The residue classes are taken in a keyed shuffle of the ranks, never of the ranks
        // themselves. A lexicographic rank is positional: contiguous blocks fixed the LEADING free
        // residues by fold (fold 0 of AEGLSVTK always put E second, fold 8 V), and residue classes
        // fixed the TRAILING ones, since rank mod 2 orders the last two free residues and rank mod 6
        // the last three (at r = 6 each fold of AEGLSVTK saw one order of them under every seed).
        // Pooled over folds that is fine, but a single fold, or an estimate per fold, drew from a
        // skewed population. A bijection keeps the shares disjoint and their sizes unchanged.
        BigInteger usable = size - BigInteger.One;
        BigInteger share = (usable - fold + foldCount - 1) / foldCount;
        byte[] key = DeriveKey(targetSequence, seed);
        BigInteger offset = ToNonNegative(key) % share;
        var shuffle = new KeyedShuffle(key, usable);

        bool anyRejectedOnlyByContext = false;
        BigInteger probes = BigInteger.Zero;

        for (BigInteger step = BigInteger.Zero; step < share; step++)
        {
            probes = step + BigInteger.One;
            BigInteger nonIdentityRank = shuffle.Apply(fold + foldCount * ((offset + step) % share));
            BigInteger index = nonIdentityRank < identity ? nonIdentityRank : nonIdentityRank + BigInteger.One;
            string candidate = DecoySequenceValidator.UnrankPermutation(targetSequence, motifs, index,
                out int[] swapped, alsoHeldInPlace);

            if (forbiddenSequences.Contains(candidate))
            {
                continue;
            }

            if (rejectInContext is not null && rejectInContext(candidate))
            {
                // Usable on its own merits, and refused only because of what it sits next to.
                anyRejectedOnlyByContext = true;
                continue;
            }

            return new EntrapmentPeptide(targetSequence, candidate, swapped, fold, size,
                Probes(probes), EntrapmentFailure.None);
        }

        // The share was walked end to end, so this is a proof rather than an abandoned search --
        // and which proof it is depends on what did the refusing. Reporting a run collision as
        // "all permutations taken" would send a caller after a different target database when the
        // answer is a different seed. The probe count is the whole share, not zero: it is the
        // evidence the walk was exhaustive, and it is on the failure paths that a reader most wants
        // to know how much was examined. The identity is never in a share, so a fold count too
        // large for the space is decided by the guard above and cannot reach here.
        EntrapmentFailure reason = anyRejectedOnlyByContext
            ? EntrapmentFailure.RunCollisionsExhaustedTheSpace
            : EntrapmentFailure.AllPermutationsTaken;

        return Failed(targetSequence, fold, size, Probes(probes), reason);
    }

    /// <summary>
    /// The probe count as a diagnostic <see cref="int"/>, saturating rather than throwing.
    /// </summary>
    /// <remarks>
    /// The walk is over a <see cref="BigInteger"/> share, which a lightly-pinned peptide can make
    /// larger than <see cref="int.MaxValue"/>. An explicit conversion is checked, so a walk that
    /// ever got that far would throw <see cref="OverflowException"/> out of a report build whose job
    /// was to classify the outcome -- trading a usable answer for a crash, over a number that is
    /// only ever read as a diagnostic. Saturating says "at least this many", which is what the
    /// column means at that magnitude anyway.
    /// </remarks>
    private static int Probes(BigInteger probes) =>
        probes >= int.MaxValue ? int.MaxValue : (int)probes;

    private static EntrapmentPeptide Failed(string targetSequence, int fold, BigInteger size,
        int probesUsed, EntrapmentFailure failure) =>
        new(targetSequence, null, null, fold, size, probesUsed, failure);

    /// <summary>
    /// The key every choice for this peptide derives from: where in a fold's share to start looking,
    /// and the shuffle the shares are taken in. Independent of the fold, so every fold of a peptide
    /// shares one shuffle and the shares stay disjoint.
    /// </summary>
    /// <remarks>
    /// SHA-256 rather than <see cref="string.GetHashCode()"/> or a home-made hash: its output is
    /// fixed by specification, so a database regenerated on another machine, another runtime or in
    /// another decade is byte-identical. String hash codes are explicitly not stable across runs.
    /// </remarks>
    private static byte[] DeriveKey(string sequence, int seed)
    {
        // Format the seed invariantly. Interpolation uses the current culture, and a negative
        // seed renders its sign as U+002D under en-US but U+2212 MINUS SIGN under sv-SE, fi-FI and
        // lt-LT -- different bytes into SHA-256, so the same request would produce a different
        // database on a differently-configured machine. The reproducibility this method exists for
        // has to survive a culture change, not only a process restart.
        string material = seed.ToString(CultureInfo.InvariantCulture) + ":" + sequence;
        return SHA256.HashData(Encoding.UTF8.GetBytes(material));
    }

    /// <summary>The bytes read as an unsigned little-endian integer.</summary>
    private static BigInteger ToNonNegative(byte[] bytes) => new(bytes, isUnsigned: true);

    /// <summary>
    /// A keyed bijection on <c>[0, domain)</c>: a Feistel network over the fewest bits that hold the
    /// domain, cycle-walked back into it.
    /// </summary>
    /// <remarks>
    /// <para>Feistel rounds are invertible whatever the round function, so this is a permutation by
    /// construction, and cycle-walking (re-applying until the value lands inside the domain) keeps it
    /// one on the domain itself. The width is the fewest bits that hold the domain, so fewer than two
    /// passes are needed on average. An odd width is split unevenly and the halves swap sizes each
    /// round, which is still invertible.</para>
    /// <para>The round function is SHA-256 over the key, the round, a block counter and the right
    /// half -- fixed by specification, like the key, so the shuffle is part of what regenerates
    /// byte-identically. An affine map <c>k -> (a k + b) mod domain</c> would be cheaper but keeps
    /// lattice structure in the low digits, which is the skew this exists to remove.</para>
    /// </remarks>
    private sealed class KeyedShuffle
    {
        private const int Rounds = 4;
        private readonly byte[] _key;
        private readonly BigInteger _domain;
        private readonly int _bits;

        internal KeyedShuffle(byte[] key, BigInteger domain)
        {
            _key = key;
            _domain = domain;
            _bits = domain <= BigInteger.One ? 0 : (int)(domain - BigInteger.One).GetBitLength();
        }

        internal BigInteger Apply(BigInteger value)
        {
            if (_bits == 0)
            {
                return value;
            }

            do
            {
                value = Encrypt(value);
            }
            while (value >= _domain);

            return value;
        }

        private BigInteger Encrypt(BigInteger value)
        {
            int replaced = _bits / 2;
            for (int round = 0; round < Rounds; round++)
            {
                int kept = _bits - replaced;
                BigInteger left = value >> kept;
                BigInteger right = value & ((BigInteger.One << kept) - BigInteger.One);
                value = (right << replaced) | (left ^ RoundFunction(round, right, replaced));
                replaced = kept;
            }
            return value;
        }

        /// <summary>A pseudo-random value of <paramref name="bits"/> bits, keyed by the round and the half.</summary>
        private BigInteger RoundFunction(int round, BigInteger right, int bits)
        {
            if (bits == 0)
            {
                return BigInteger.Zero;
            }

            int rightLength = right.GetByteCount(isUnsigned: true);
            byte[] input = new byte[_key.Length + 2 + rightLength];
            _key.CopyTo(input, 0);
            input[_key.Length] = (byte)round;
            right.TryWriteBytes(input.AsSpan(_key.Length + 2), out _, isUnsigned: true);

            int blocks = (bits + 255) / 256;
            byte[] stream = new byte[blocks * 32];
            for (int block = 0; block < blocks; block++)
            {
                input[_key.Length + 1] = (byte)block;
                SHA256.HashData(input, stream.AsSpan(block * 32, 32));
            }

            return new BigInteger(stream, isUnsigned: true) & ((BigInteger.One << bits) - BigInteger.One);
        }
    }
}
