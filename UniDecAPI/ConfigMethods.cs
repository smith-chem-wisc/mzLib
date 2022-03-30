using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using MassSpectrometry; 

namespace UniDecAPI
{
	public partial class UniDecAPIMethods
	{
		public static class ConfigMethods
		{
			[DllImport("TestDLL.dll", EntryPoint = "ReturnModifiedConfigToCS")]
			private static extern Config _ReturnModifiedConfigToCS(Config config); 
			public static Config ReturnModifiedConfigToCS(Config config)
			{
				return _ReturnModifiedConfigToCS(config); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "ModifyConfigToDefault")]
			private static unsafe extern Config _ModifyConfigToDefault(Config* config); 
			public static unsafe Config ModifyConfigToDefault(Config config)
			{
				Config* ptrToConfig = &config; 
				return _ModifyConfigToDefault(ptrToConfig); 
			}

			[DllImport("TestDLL.dll", EntryPoint = "PostImport")]
			private static unsafe extern void _PostImport(ref Config config);
			public static void PostImport(ref Config config)
			{
				unsafe
				{
					_PostImport(ref config);
				}
			}

			public static void SetupConfig(ref Config config, MsDataScan scan)
			{
				config.lengthmz = scan.MassSpectrum.XArray.Length;
				config.filterwidth = 10;
				config.minmz = (float)scan.MassSpectrum.XArray.Min();
				config.maxmz = (float)scan.MassSpectrum.XArray.Max();
			}
			public static void CreateAndSetupConfig(MsDataScan scan, out Config config)
			{
				config = new(); 
				config = ModifyConfigToDefault(config);
				SetupConfig(ref config, scan);
				PostImport(ref config);
			}
		}
	}
}
