using System;
using System.Numerics;

namespace StatisticalModels
{
    /// <summary>
    /// Feature rows of doubles or floats. A float widens to a double exactly, so a model reads a float row as the same values
    /// held in a double row, and scores it bit for bit alike.
    /// </summary>
    internal static class FeatureRows
    {
        /// <summary>The row itself when it is a double row; otherwise its values widened into a new one.</summary>
        public static double[] Widened<T>(T[] row) where T : unmanaged, IBinaryFloatingPointIeee754<T>
        {
            if (row is double[] doubles)
                return doubles;
            var widened = new double[row.Length];
            for (int j = 0; j < row.Length; j++)
                widened[j] = double.CreateTruncating(row[j]);
            return widened;
        }

        /// <summary>Whether every value is finite.</summary>
        public static bool AllFinite<T>(T[] row) where T : unmanaged, IBinaryFloatingPointIeee754<T>
        {
            foreach (T value in row)
                if (!T.IsFinite(value))
                    return false;
            return true;
        }
    }
}
