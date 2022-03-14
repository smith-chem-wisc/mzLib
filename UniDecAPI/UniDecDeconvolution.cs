using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry; 

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
				config.lengthmz = xarray.Length; 
				Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs(); 
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					try
					{
						inp.dataMZ = ptrXarray;
						inp.dataInt = ptrYarray;

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
		public static int TestPerformUniDecDecon(this MsDataScan scan)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			int result = 0; 
			unsafe
			{
				Config config = new();
				config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
				// Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
						inp.dataMZ = ptrXarray;
						inp.dataInt = ptrYarray;
						result = UniDecAPIMethods.DeconMethods.TestingMainDeconvolution();
				}
			}
			return result; 
		}
	}
}
