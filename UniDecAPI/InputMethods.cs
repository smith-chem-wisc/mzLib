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
		public static class InputMethods
		{
			[DllImport("TestDLL.dll", EntryPoint = "SetupInputs")]
			private static extern Input _SetupInputs(); 
			public static Input SetupInputs()
			{
				return _SetupInputs(); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "FreeInputs")]
			private static extern void _FreeInputs(Input inp); 
			public static void FreeInputs(Input inp)
			{
				_FreeInputs(inp); 
			}
		}
	}
}
