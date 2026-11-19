using System;

namespace SpectralAveraging;

/// <summary>
///     Working arrays for rejecting outliers bin after bin without allocating for every bin. Not thread-safe: each thread
///     needs its own.
/// </summary>
internal sealed class OutlierRejectionBuffers
{
    private double[] _values = new double[16];
    private double[][] _exactLength = new double[17][];

    /// <summary>
    ///     An array of at least <paramref name="minimumLength"/> elements; its contents are not kept between calls.
    /// </summary>
    internal double[] Values(int minimumLength)
    {
        if (_values.Length < minimumLength)
            _values = new double[Math.Max(minimumLength, 2 * _values.Length)];
        return _values;
    }

    /// <summary>
    ///     An array of exactly <paramref name="length"/> elements, for MathNet methods that work on a whole array; its
    ///     contents are not kept between calls.
    /// </summary>
    internal double[] ExactLength(int length)
    {
        if (length >= _exactLength.Length)
            Array.Resize(ref _exactLength, Math.Max(length + 1, 2 * _exactLength.Length));
        return _exactLength[length] ??= new double[length];
    }
}
