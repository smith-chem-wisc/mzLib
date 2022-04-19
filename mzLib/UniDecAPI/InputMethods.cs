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
			[DllImport("TestDLL.dll", EntryPoint = "ReadInputs")]
			private static extern void _ReadInputs(ref InputUnsafe inp, ref Config config);
			public static void ReadInputs(InputUnsafe inp, Config config)
			{
				_ReadInputs(ref inp, ref config); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "ReadInputsByValue")]
			private static extern InputUnsafe _ReadInputsByValue(InputUnsafe inp, ref Config config); 
			public static InputUnsafe ReadInputsByValue(InputUnsafe inp, Config config)
			{
				return _ReadInputsByValue(inp, ref config); 
			}
		}
	}
}
