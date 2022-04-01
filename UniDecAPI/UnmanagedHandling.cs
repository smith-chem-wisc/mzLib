using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public static class UnmanagedHandling
	{
		// takes a managed array and converts it to an IntPtr. 
		// You need to be sure that you free any IntPtr you create with this method. 
		// Immediately castable to whatever you want your pointer to be. 
		public static IntPtr AllocateToUnmanaged(float[] managedArray)
		{
			IntPtr ptr = Marshal.AllocHGlobal(Marshal.SizeOf(managedArray[0]) * managedArray.Length);
			Marshal.Copy(managedArray, 0, ptr, managedArray.Length);
			return ptr; 
		}
		public static IntPtr AllocateToUnmanaged(int size)
		{
			IntPtr ptr = Marshal.AllocHGlobal(size);
			return ptr; 
		}
	}
}
