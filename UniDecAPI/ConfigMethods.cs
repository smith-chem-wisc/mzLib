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
		public static class ConfigMethods
		{
			[DllImport("UniDecMinimal", EntryPoint = "SetDefaultConfig")]
			private static extern Config _SetDefaultConfig();
			public static Config SetDefaultConfig()
			{
				return SetDefaultConfig();
			}
		}
	}
}
