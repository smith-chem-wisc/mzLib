using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UniDecAPI
{
	public static class ArrayConversionExtensions
	{
		public static float[] ConvertDoubleArrayToFloat(this double[] array)
		{
			float[] floatArray = new float[array.Length];
			for (int i = 0; i < array.Length; i++)
			{
				floatArray[i] = (float)array[i];
			}
			return floatArray;
		}
	}
}
