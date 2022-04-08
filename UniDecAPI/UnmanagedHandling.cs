using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public class UnmanagedHandling : IDisposable
	{
		// TODO: Implement IDisposable.

		public List<IntPtr> ListPtrs { get; private set; }
		public UnmanagedHandling()
		{
			ListPtrs = new List<IntPtr>(); 
		}

		// takes a managed array and converts it to an IntPtr. 
		// You need to be sure that you free any IntPtr you create with this method. 
		// Immediately castable to whatever you want your pointer to be. 
		public IntPtr AllocateToUnmanaged(float[] managedArray)
		{
			IntPtr ptr = Marshal.AllocHGlobal(Marshal.SizeOf(managedArray[0]) * managedArray.Length);
			Marshal.Copy(managedArray, 0, ptr, managedArray.Length);
			ListPtrs.Add(ptr); 
			return ptr; 
		}
		public IntPtr AllocateToUnmanaged(int size, Type type)
		{
			IntPtr ptr = Marshal.AllocHGlobal(size * Marshal.SizeOf(type));
			ListPtrs.Add(ptr); 
			return ptr; 
		}
		public void Dispose()
		{
			Dispose(true);
			GC.SuppressFinalize(this);
		}
		protected virtual void Dispose(bool disposing)
		{
			if (disposing)
			{

			}
			foreach(IntPtr ptrList in ListPtrs)
			{
				Marshal.FreeHGlobal(ptrList);
			}
		}
	}
}
