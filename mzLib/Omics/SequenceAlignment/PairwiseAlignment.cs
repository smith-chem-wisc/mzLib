using System;
using System.Collections.Generic;
using System.Text;

namespace Omics.SequenceAlignment;

/// <summary>
/// One alignment of two sequences, from <see cref="PairwiseAligner.Align"/>: which residue of each
/// corresponds to which residue of the other, and which residues face a gap.
/// </summary>
/// <remarks>
/// Positions are 1-based, as everywhere else in mzLib (<c>OneBasedPossibleLocalizedModifications</c>,
/// <c>OneBasedBeginPosition</c>). A residue facing a gap has no counterpart, and the Try methods
/// return false for it rather than the nearest residue.
/// </remarks>
public sealed class PairwiseAlignment
{
    private readonly int[] _sourceToTarget;
    private readonly int[] _targetToSource;

    internal PairwiseAlignment(string source, string target, IReadOnlyList<(int Source, int Target)> columns,
        int score, string alignerId)
    {
        Source = source;
        Target = target;
        Score = score;
        AlignerId = alignerId;
        _sourceToTarget = new int[source.Length + 1];
        _targetToSource = new int[target.Length + 1];

        var alignedSource = new StringBuilder(columns.Count);
        var alignedTarget = new StringBuilder(columns.Count);
        foreach (var (s, t) in columns)
        {
            alignedSource.Append(s == 0 ? '-' : source[s - 1]);
            alignedTarget.Append(t == 0 ? '-' : target[t - 1]);
            if (s != 0 && t != 0)
            {
                _sourceToTarget[s] = t;
                _targetToSource[t] = s;
                if (char.ToUpperInvariant(source[s - 1]) == char.ToUpperInvariant(target[t - 1]))
                {
                    Identical++;
                }
                else
                {
                    Substituted++;
                }
            }
            else if (s != 0)
            {
                SourceResiduesFacingGap++;
            }
            else
            {
                TargetResiduesFacingGap++;
            }
        }
        AlignedSource = alignedSource.ToString();
        AlignedTarget = alignedTarget.ToString();
    }

    /// <summary>The first sequence, as given.</summary>
    public string Source { get; }

    /// <summary>The second sequence, as given.</summary>
    public string Target { get; }

    /// <summary>The source with <c>-</c> where it faces a target residue across a gap; as long as <see cref="AlignedTarget"/>.</summary>
    public string AlignedSource { get; }

    /// <summary>The target with <c>-</c> where it faces a source residue across a gap.</summary>
    public string AlignedTarget { get; }

    /// <summary>The alignment's score under the aligner that made it.</summary>
    public int Score { get; }

    /// <summary>The aligner and its parameters (<see cref="PairwiseAligner.Id"/>).</summary>
    public string AlignerId { get; }

    /// <summary>Paired residues that are the same amino acid (compared case-insensitively).</summary>
    public int Identical { get; }

    /// <summary>Paired residues that differ.</summary>
    public int Substituted { get; }

    /// <summary>Source residues with no target counterpart.</summary>
    public int SourceResiduesFacingGap { get; }

    /// <summary>Target residues with no source counterpart.</summary>
    public int TargetResiduesFacingGap { get; }

    /// <summary>The target residue paired with a source residue.</summary>
    /// <param name="oneBasedSourcePosition">A position in <see cref="Source"/>, 1 to its length.</param>
    /// <param name="oneBasedTargetPosition">The paired target position, or 0 when the residue faces a gap.</param>
    /// <returns>False when the source residue faces a gap.</returns>
    /// <exception cref="ArgumentOutOfRangeException">The position is outside the source.</exception>
    public bool TryMapSourceToTarget(int oneBasedSourcePosition, out int oneBasedTargetPosition)
        => TryMap(_sourceToTarget, oneBasedSourcePosition, out oneBasedTargetPosition);

    /// <summary>The source residue paired with a target residue.</summary>
    /// <param name="oneBasedTargetPosition">A position in <see cref="Target"/>, 1 to its length.</param>
    /// <param name="oneBasedSourcePosition">The paired source position, or 0 when the residue faces a gap.</param>
    /// <returns>False when the target residue faces a gap.</returns>
    /// <exception cref="ArgumentOutOfRangeException">The position is outside the target.</exception>
    public bool TryMapTargetToSource(int oneBasedTargetPosition, out int oneBasedSourcePosition)
        => TryMap(_targetToSource, oneBasedTargetPosition, out oneBasedSourcePosition);

    private static bool TryMap(int[] map, int position, out int mapped)
    {
        if (position < 1 || position >= map.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(position), position, $"Expected 1 to {map.Length - 1}.");
        }
        mapped = map[position];
        return mapped != 0;
    }
}
