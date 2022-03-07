using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace TestUniDec
{
	[StructLayout(LayoutKind.Sequential)]
	public struct TestStruct
	{
		public int val1;
	}
	[StructLayout(LayoutKind.Sequential)]
	public struct TestStruct2
	{
		public int val1;
		public IntPtr array1; 
	}

	// Shows how to convert a struct to a class using a parameterized constructor taking unmanaged code matching struct
	public class Struct2
	{
		public int Val1 { get; set; }
		public int[] Array1 { get; set; }
		public Struct2(TestStruct2 ts2, int array1Length)
		{
			// Use the constructor to do type marshalling for each pointer that the unmanaged 
			// code returns. 
			Val1 = ts2.val1;
			// This is an idiom that you'll need to know. 
			int[] _array1 = new int[array1Length];
			Marshal.Copy(ts2.array1, _array1, 0, array1Length);
			Array1 = _array1; 
		}
		public Struct2()
		{

		}
	}

	public unsafe struct TestStruct2Unsafe
	{
		public int val1;
		public int* array1; 
	}
}
