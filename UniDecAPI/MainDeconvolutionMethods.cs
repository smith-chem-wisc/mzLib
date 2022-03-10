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
		public static class MainDeconvolutionMethods
		{
			[DllImport("UniDecMinimal", EntryPoint = "MainDeconvolution")]
			private static extern Decon _MainDeconvolution(Config config, Input inp, int silent, int verbose);
			public static void MainDeconvolution(Config config, Input inp, int silent, int verbose, out Decon decon)
			{
				try
				{ 
					decon = _MainDeconvolution(config, inp, silent, verbose);
				}
				finally
				{
					InputMethods.FreeInputs(inp); 
				}
			} 
		}
	}
}
