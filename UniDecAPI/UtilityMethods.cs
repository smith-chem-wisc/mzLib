using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 
namespace UniDecAPI
{
	public partial class UniDecAPIMethods
	{
		public class UtilityMethods
		{
			public static unsafe void KillB_CS(ref float[] I, bool[] B, float intthresh, int lengthmz, int numz, int isolength, int[] isotopeops, float[] isotopeval)
			{
				int i, j, k; 
				if(isolength == 0)
				{
					for(i = 0; i < lengthmz; i++)
					{
						for(j = 0; j < numz; j++)
						{

						}
					}
				}
				else
				{
					float cutoff = 0.5F; 
					for(i = 0; i < lengthmz; i++)
					{
						for(j = 0; j < numz; j++)
						{
							float max = 0F; 
							for(k = 0; k < isolength; k++)
							{
								float val = isotopeval[index3D(numz, isolength, i, j, k)]; 
								if (val > max) { max = val;  }
							}
							for(k = 0; k < isolength; k++)
							{
								float val = isotopeval[index3D(numz, isolength, i, j, k)]; 
								if(val > cutoff * max)
								{
									int pos = isotopeops[index3D(numz, isolength, i, j, k)]; 
									if(I[pos] <= intthresh) { B[index2D(numz, i, j)] = false;  }
								}
							}
						}
					}
				}
			}
			public static int index3D(int ncols, int nrows, int r, int c, int d)
			{
				return r * ncols * nrows + c * nrows + d;
			}
			public static int index2D(int ncols, int r, int c)
			{
				return r * ncols + c;
			}
			// Creating an array from a pointer is not very fast, so you need to 
			// use sparingly. 
			public static unsafe float[] PtrToArray(float* pointer, int length)
			{
				float[] result = new float[length]; 
				for(int i = 0; i < length; i++)
				{
					result[i] = pointer[i]; 
				}
				return result; 
			}
			public static unsafe char[] PtrToArray(char* pointer, int length)
			{
				char[] result = new char[length]; 
				for(int i = 0; i < length; i++)
				{
					result[i] = pointer[i]; 
				}
				return result;
			}
			public static unsafe byte[] PtrToArray(byte* pointer, int length)
			{
				byte[] result = new byte[length];
				for (int i = 0; i < length; i++)
				{
					result[i] = pointer[i]; 					
				}
				return result;
			}
			public static void ConvertASCIIBytesFromCToByteInCS(ref byte[] byteArray)
			{
				for(int i = 0; i < byteArray.Length; i++)
				{
					// 48 is 0 in ascii and 57 is 9 in ascii
					if (byteArray[i] >= 48 & byteArray[i] <= 57)
					{
						byteArray[i] -= (byte)48;
					}
				}
			}
			public static void ConvertASCIIBytesFromCToByteInCS(ref byte byteVal)
			{
					// 48 is 0 in ascii and 57 is 9 in ascii
					if (byteVal >= 48 & byteVal <= 57)
					{
						byteVal -= (byte)48;
					}
			}
			[DllImport("TestDLL.dll", EntryPoint = "SetLimits")]
			private static extern unsafe void _SetLimits(Config config, InputUnsafe* inp); 
			public static unsafe void SetLimits(Config config, InputUnsafe inp)
			{
				_SetLimits(config, &inp); 
			}

		}
	}
}
