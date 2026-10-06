using System;
using System.Collections.Generic;
using System.Globalization;

namespace Omics.SequenceAlignment;

/// <summary>
/// Global pairwise alignment of two residue sequences with affine gaps (Gotoh's form of
/// Needleman-Wunsch), giving the residue-to-residue correspondence between them.
/// </summary>
/// <remarks>
/// <para>
/// A gap of length <c>k</c> costs <c>GapOpen + k * GapExtend</c>, which is BLAST's convention:
/// with the defaults (BLOSUM62, 11, 1) a one-residue gap costs 12.
/// </para>
/// <para>
/// With <see cref="FreeEndGaps"/> (the default) a gap at either end of either sequence costs
/// nothing, so one sequence may extend past the other at either terminus without penalty: an
/// entry with a signal peptide the other lacks, or a longer C-terminus. EMBOSS needle scores end
/// gaps the same way by default. With it off, every gap is charged.
/// </para>
/// <para>
/// The result is one optimal alignment. Where several alignments score the same, the traceback
/// prefers, at each step, a residue pair over a gap in the target over a gap in the source, so the
/// same inputs always give the same alignment. Two identical sequences are aligned residue for
/// residue without the dynamic programme when <see cref="SubstitutionMatrix.IdenticalSequencesAlignToThemselves"/>
/// shows that alignment is optimal, which it is for BLOSUM62 and sequences free of <c>X</c> and <c>*</c>.
/// </para>
/// <para>
/// Time is proportional to the product of the two lengths, and so is memory: one byte per cell
/// for the traceback. <see cref="MaxCells"/> bounds it, and a larger pair is refused rather than
/// allowed to exhaust memory.
/// </para>
/// </remarks>
public sealed class PairwiseAligner
{
    private const int NegativeInfinity = int.MinValue / 4;

    // Traceback bits, per cell: which state each of the three states came from.
    private const int FromMatch = 0, FromSourceGap = 1, FromTargetGap = 2;

    /// <summary>Creates an aligner.</summary>
    /// <param name="matrix">The substitution matrix; <see cref="SubstitutionMatrix.Blosum62"/> when null.</param>
    /// <param name="gapOpen">The cost of opening a gap, charged once per gap. Not negative.</param>
    /// <param name="gapExtend">The cost of each residue in a gap, including the first. Greater than 0.</param>
    /// <param name="freeEndGaps">Whether a gap at either end of either sequence is free.</param>
    /// <param name="maxCells">The largest <c>(source length + 1) x (target length + 1)</c> aligned with the dynamic programme.</param>
    public PairwiseAligner(SubstitutionMatrix? matrix = null, int gapOpen = 11, int gapExtend = 1,
        bool freeEndGaps = true, long maxCells = 200_000_000)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(gapOpen);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(gapExtend);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxCells);
        Matrix = matrix ?? SubstitutionMatrix.Blosum62;
        GapOpen = gapOpen;
        GapExtend = gapExtend;
        FreeEndGaps = freeEndGaps;
        MaxCells = maxCells;
    }

    /// <summary>The substitution matrix.</summary>
    public SubstitutionMatrix Matrix { get; }

    /// <summary>The cost of opening a gap.</summary>
    public int GapOpen { get; }

    /// <summary>The cost of each residue in a gap.</summary>
    public int GapExtend { get; }

    /// <summary>Whether gaps at the ends of either sequence are free.</summary>
    public bool FreeEndGaps { get; }

    /// <summary>The largest dynamic-programming table, in cells, that <see cref="Align"/> will build.</summary>
    public long MaxCells { get; }

    /// <summary>
    /// Names the method and every parameter, for example
    /// <c>global-affine;BLOSUM62;open=11;extend=1;end-gaps=free</c>, so a stored correspondence can
    /// say exactly how it was aligned.
    /// </summary>
    public string Id => string.Create(CultureInfo.InvariantCulture,
        $"global-affine;{Matrix.Name};open={GapOpen};extend={GapExtend};end-gaps={(FreeEndGaps ? "free" : "charged")}");

    /// <summary>Aligns <paramref name="source"/> with <paramref name="target"/>.</summary>
    /// <param name="source">The first sequence, as single-letter residue codes.</param>
    /// <param name="target">The second sequence.</param>
    /// <exception cref="ArgumentNullException">Either sequence is null.</exception>
    /// <exception cref="ArgumentException">The table would exceed <see cref="MaxCells"/>.</exception>
    public PairwiseAlignment Align(string source, string target)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentNullException.ThrowIfNull(target);

        if (string.Equals(source, target, StringComparison.OrdinalIgnoreCase)
            && Matrix.IdenticalSequencesAlignToThemselves && !ContainsUnknownOrStop(source))
        {
            var pairs = new List<(int, int)>(source.Length);
            int score = 0;
            for (int k = 1; k <= source.Length; k++)
            {
                pairs.Add((k, k));
                score += Matrix.Score(source[k - 1], target[k - 1]);
            }
            return new PairwiseAlignment(source, target, pairs, score, Id);
        }

        int n = source.Length, m = target.Length;
        long cells = (long)(n + 1) * (m + 1);
        if (cells > MaxCells)
        {
            throw new ArgumentException(
                $"Aligning lengths {n} and {m} needs {cells:N0} cells, more than MaxCells ({MaxCells:N0}).");
        }
        return Gotoh(source, target, n, m);
    }

    private bool ContainsUnknownOrStop(string sequence)
    {
        foreach (char c in sequence)
        {
            if (Matrix.IsUnknownOrStop(c))
            {
                return true;
            }
        }
        return false;
    }

    private PairwiseAlignment Gotoh(string source, string target, int n, int m)
    {
        int[] s = new int[n], t = new int[m];
        for (int i = 0; i < n; i++) s[i] = Matrix.IndexOf(source[i]);
        for (int j = 0; j < m; j++) t[j] = Matrix.IndexOf(target[j]);

        int open = GapOpen + GapExtend, extend = GapExtend;
        int width = m + 1;
        // Bits 0-1: where Match came from; 2-3: SourceGap; 4-5: TargetGap.
        var trace = new byte[(n + 1) * width];

        // Match: source[i] with target[j]. SourceGap: source[i] with a gap. TargetGap: target[j] with a gap.
        var prevM = new int[width]; var prevS = new int[width]; var prevT = new int[width];
        var curM = new int[width]; var curS = new int[width]; var curT = new int[width];

        prevM[0] = 0; prevS[0] = NegativeInfinity; prevT[0] = NegativeInfinity;
        for (int j = 1; j <= m; j++)
        {
            prevM[j] = NegativeInfinity;
            prevS[j] = NegativeInfinity;
            prevT[j] = FreeEndGaps ? 0 : -(GapOpen + j * GapExtend);
            trace[j] = FromTargetGap << 4;
        }

        // Free trailing gaps: the best end on the last column (i, m) for every i.
        var lastColumn = new (int Score, int State)[n + 1];
        lastColumn[0] = Best(prevM[m], prevS[m], prevT[m]);

        for (int i = 1; i <= n; i++)
        {
            curM[0] = NegativeInfinity;
            curS[0] = FreeEndGaps ? 0 : -(GapOpen + i * GapExtend);
            curT[0] = NegativeInfinity;
            trace[i * width] = FromSourceGap << 2;

            int si = s[i - 1];
            for (int j = 1; j <= m; j++)
            {
                // Match from the diagonal.
                var (diag, diagFrom) = Best(prevM[j - 1], prevS[j - 1], prevT[j - 1]);
                curM[j] = diag == NegativeInfinity ? NegativeInfinity : diag + Matrix.ScoreByIndex(si, t[j - 1]);

                // Source residue against a gap: from the cell above.
                var (up, upFrom) = Best(prevM[j] - open, prevS[j] - extend, prevT[j] - open);
                curS[j] = up;

                // Target residue against a gap: from the cell to the left.
                var (left, leftFrom) = Best(curM[j - 1] - open, curS[j - 1] - open, curT[j - 1] - extend);
                curT[j] = left;

                trace[i * width + j] = (byte)(diagFrom | (upFrom << 2) | (leftFrom << 4));
            }
            lastColumn[i] = Best(curM[m], curS[m], curT[m]);

            (prevM, curM) = (curM, prevM);
            (prevS, curS) = (curS, prevS);
            (prevT, curT) = (curT, prevT);
        }

        // prev* now hold the last row (i = n). Choose where the alignment ends.
        int endI = n, endJ = m;
        var (bestScore, bestState) = Best(prevM[m], prevS[m], prevT[m]);
        if (FreeEndGaps)
        {
            // Prefer the full corner on a tie, then the longest prefix of each sequence.
            for (int j = m - 1; j >= 0; j--)
            {
                var (score, state) = Best(prevM[j], prevS[j], prevT[j]);
                if (score > bestScore) { bestScore = score; bestState = state; endI = n; endJ = j; }
            }
            for (int i = n - 1; i >= 0; i--)
            {
                if (lastColumn[i].Score > bestScore) { (bestScore, bestState) = lastColumn[i]; endI = i; endJ = m; }
            }
        }

        var reversed = new List<(int, int)>(n + m);
        for (int i = n; i > endI; i--) reversed.Add((i, 0));
        for (int j = m; j > endJ; j--) reversed.Add((0, j));

        int ci = endI, cj = endJ, stateNow = bestState;
        while (ci > 0 || cj > 0)
        {
            int bits = trace[ci * width + cj];
            if (ci == 0) stateNow = FromTargetGap;
            else if (cj == 0) stateNow = FromSourceGap;

            switch (stateNow)
            {
                case FromMatch:
                    reversed.Add((ci, cj));
                    stateNow = bits & 3;
                    ci--; cj--;
                    break;
                case FromSourceGap:
                    reversed.Add((ci, 0));
                    stateNow = (bits >> 2) & 3;
                    ci--;
                    break;
                default:
                    reversed.Add((0, cj));
                    stateNow = (bits >> 4) & 3;
                    cj--;
                    break;
            }
        }
        reversed.Reverse();
        return new PairwiseAlignment(source, target, reversed, bestScore, Id);
    }

    // Highest score; ties go to Match, then SourceGap, then TargetGap.
    private static (int Score, int State) Best(int match, int sourceGap, int targetGap)
    {
        match = Math.Max(match, NegativeInfinity);
        sourceGap = Math.Max(sourceGap, NegativeInfinity);
        targetGap = Math.Max(targetGap, NegativeInfinity);
        if (match >= sourceGap && match >= targetGap) return (match, FromMatch);
        return sourceGap >= targetGap ? (sourceGap, FromSourceGap) : (targetGap, FromTargetGap);
    }
}
