using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static class UniDecDeconvolution
	{
		public static Decon PerformUniDecDeconvolution(this MsDataScan scan, int silent, int verbose)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			unsafe
			{
				Config config = new();
				config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
				UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);

				Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				char[] test = { (char)0, (char)0, (char)1 };
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					fixed (char* ptrToTest = &test[0])
					{
						try
						{

							inp.dataMZ = ptrXarray;
							inp.dataInt = ptrYarray;
							inp.barr = ptrToTest;

							decon = UniDecAPIMethods.DeconMethods.MainDeconvolution(config, inp, silent, verbose);
						}
						finally
						{
							UniDecAPIMethods.InputMethods.FreeInputs(inp);
						}
						return decon;
					}
				}
			}

		}
		public static Decon TestPerformUniDecDecon(this MsDataScan scan)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			int result = 0;
			unsafe
			{
				Config config = new();
				config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
				UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);

				Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				char[] test = { (char)0, (char)0, (char)1 };
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					fixed (char* ptrToTest = &test[0])
					{
						inp.dataMZ = ptrXarray;
						inp.dataInt = ptrYarray;
						inp.barr = ptrToTest;
						decon = UniDecAPIMethods.DeconMethods.TestingMainDeconvolution(config, inp, 0, 1);
						
						return decon;
					}
				}
			}
		}
	}
}
