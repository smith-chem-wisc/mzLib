using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static class ArrayIndexing
	{
		public static int Index2D(int ncols, int r, int c)
		{
			return r * ncols + c;
		}
		public static int Index3D(int ncols, int nrows, int r, int c, int d)
		{
			return r * ncols * nrows + c * nrows + d;
		}
		public static int Indexmod(int length, int r, int c)
		{
			int a = c - r;
			int b = length;
			int result = a % b;
			return result < 0 ? result + b : result;
		}
		public static void ApplyCutoff1D(ref float[] array, float cutoff, int lengthmz)
		{
			for (int i = 0; i < lengthmz; i++)
			{
				if (array[i] < cutoff)
				{
					array[i] = 0;
				}
			}
		}
	}
}
