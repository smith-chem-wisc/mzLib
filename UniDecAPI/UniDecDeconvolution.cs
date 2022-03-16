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
		public static Decon PerformUniDecDeconvolution(this MsDataScan scan)
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
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{
					// need to fully setup the config object. 
						try
						{

							inp.dataMZ = ptrXarray;
							inp.dataInt = ptrYarray;

							decon = UniDecAPIMethods.DeconMethods.RunUnidec(inp,config);
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
}
