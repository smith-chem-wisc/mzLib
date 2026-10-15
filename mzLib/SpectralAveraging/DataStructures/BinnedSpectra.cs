using System;
using System.Buffers;

namespace SpectralAveraging;

/// <summary>
///     The peaks of a group of spectra sorted into m/z bins, all of them in one array. Bins are numbered in the order
///     they are first met, spectrum by spectrum and peak by peak - the order a <c>Dictionary</c> keyed on the bin index
///     enumerates them in. Each bin's peaks are contiguous: its own peaks in spectrum order, then a zero-intensity peak
///     for every spectrum that has none in the bin. The arrays are rented from <see cref="ArrayPool{T}"/>; dispose to
///     return them.
/// </summary>
internal sealed class BinnedSpectra : IDisposable
{
    private BinnedPeak[] _peaks;
    private int[] _start;
    private int[] _length;

    internal BinnedSpectra(BinnedPeak[] peaks, int[] start, int[] length, int count)
    {
        _peaks = peaks;
        _start = start;
        _length = length;
        Count = count;
    }

    /// <summary>
    ///     Number of bins.
    /// </summary>
    internal int Count { get; }

    /// <summary>
    ///     The peaks of the <paramref name="bin"/>-th bin to be met.
    /// </summary>
    internal Span<BinnedPeak> this[int bin] => _peaks.AsSpan(_start[bin], _length[bin]);

    public void Dispose()
    {
        if (_peaks is null)
            return;
        ArrayPool<BinnedPeak>.Shared.Return(_peaks);
        ArrayPool<int>.Shared.Return(_start);
        ArrayPool<int>.Shared.Return(_length);
        _peaks = null;
        _start = null;
        _length = null;
    }
}
