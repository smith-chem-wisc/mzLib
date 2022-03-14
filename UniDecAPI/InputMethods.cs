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
		public static class InputMethods
		{
			[DllImport("TestDLL.dll", EntryPoint = "SetupInputs")]
			private static extern InputUnsafe _SetupInputs(); 
			public static InputUnsafe SetupInputs()
			{
				return _SetupInputs(); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "FreeInputs")]
			private static extern void _FreeInputs(InputUnsafe inp); 
			public static void FreeInputs(InputUnsafe inp)
			{
				_FreeInputs(inp); 
			}
		}
	}
}
