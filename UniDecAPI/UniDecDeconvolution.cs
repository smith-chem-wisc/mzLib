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
		public static Decon PerformUniDecDeconvolution(this MsDataScan scan, Config config)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			unsafe
			{
				UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);

				Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					// need to fully setup the config object. 
					try
					{

						inp.dataMZ = ptrXarray;
						inp.dataInt = ptrYarray;
						inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config); 

						decon = UniDecAPIMethods.DeconMethods.RunUnidec(inp,config);
					}
					finally
					{

					}
					return decon;
				}
			}

		}
		public static Decon PerformUniDecDeconForTestingPurposes(this MsDataScan scan, Config config)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			unsafe
			{
				UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);

				Decon decon = new();
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					// need to fully setup the config object. 
					try
					{

						inp.dataMZ = ptrXarray;
						inp.dataInt = ptrYarray;
						inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config); 

						decon = AlgorithmTesting.RunUniDecWithTestMainDeconAlgo(inp, config);
					}
					finally
					{
						//UniDecAPIMethods.InputMethods.FreeInputs(inp);
					}
					return decon;
				}
			}
		}
		public static unsafe void Deconvolute(this MsDataScan scan, Config config)
		{
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);
			InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
			fixed (float* ptrXarry = &xarray[0], ptrYarray = &yarray[0])
			{
				inp.dataInt = ptrYarray;
				inp.dataMZ = ptrXarry;
				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config); 
			}
		}
	}
	public static class UniDecDeconvolutionUtilities
	{

	}
}
