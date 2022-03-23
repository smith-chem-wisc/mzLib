using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public static class AlgorithmTesting
	{
		[DllImport("TestDLL.dll", EntryPoint = "MemoryAllocationOfBarr")]
		private static extern int _MemoryAllocationOfBarr();
		public static int MemoryAllocationOfBarr()
		{
			return _MemoryAllocationOfBarr(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "AllocateMemoryToPointersThenFree")]
		private static extern int _AllocateMemoryToPointersThenFree(); 
		public static int AllocateMemoryToPointersThenFree()
		{
			return _AllocateMemoryToPointersThenFree(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "UseMemcpyInC")]
		private static extern char _UseMemcpyInC(); 
		public static char UseMemcpyInC()
		{
			return _UseMemcpyInC(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "UseMemcpyWithInpAndConfigObjects")]
		private static extern char _UseMemcpyWithInpAndConfigObjects(); 
		public static char UseMemcpyWithInpAndConfigObjects()
		{
			return _UseMemcpyWithInpAndConfigObjects(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "MemoryObjectAllocationToHeap")]
		private static extern int _MemoryObjectAllocationToHeap(Config config, InputUnsafe inp);
		public static int MemoryObjectAllocationToHeap(Config config, InputUnsafe inp)
		{
			return _MemoryObjectAllocationToHeap(config, inp); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "TestSetStartEnds")]
		private static extern int _TestSetStartEnds(InputUnsafe inp, Config config);
		public static int TestSetStartEnds(InputUnsafe inp, Config config)
		{
			return _TestSetStartEnds(inp, config); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "MemoryObjectAllocationToHeapConfigPtr")]
		private static extern int _MemoryObjectAllocationToHeapConfigPtr(ref Config config, InputUnsafe inp); 
		public static int MemoryObjectAllocationToHeapConfigPtr(Config config, InputUnsafe inp)
		{
			return _MemoryObjectAllocationToHeapConfigPtr(ref config, inp); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "TestFreeDecon")]
		private static extern int _TestFreeDecon(); 
		public static int TestFreeDecon()
		{
			return _TestFreeDecon(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "TestSetupAndAllocateMemoryToDecon")]
		private static extern int _TestSetupAndAllocateMemoryToDecon();
		public static int TestSetupAndAllocateMemoryToDecon()
		{
			return _TestSetupAndAllocateMemoryToDecon(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "TestSetupAndReturnDecon")]
		private static extern Decon _TestSetupAndReturnDecon();
		public static Decon TestSetupAndReturnDecon()
		{
			return _TestSetupAndReturnDecon();
		}

		[DllImport("TestDLL.dll", EntryPoint = "MainDeconWithMinimalControlFlow")]
		private static extern Decon _MainDeconWithMinimalControlFlow(Config config, InputUnsafe inp); 
		public static Decon MainDeconWithMinimalControlFlow(Config config, InputUnsafe inp)
		{
			return _MainDeconWithMinimalControlFlow(config, inp); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "RunUniDecWithTestMainDeconAlgo")]
		private static extern Decon _RunUniDecWithTestMainDeconAlgo(InputUnsafe inp, Config config); 
		public static Decon RunUniDecWithTestMainDeconAlgo(InputUnsafe inp, Config config)
		{
			return _RunUniDecWithTestMainDeconAlgo(inp, config); 	
		}
		[DllImport("TestDLL.dll", EntryPoint = "TestingKillBFunction")]
		private static extern unsafe void _TestKillBFunction(float* I, byte* B, float intthresh, int lengthmz, int numz, int isolength, int* isotopepos, float* isotopeval);
		public static unsafe void TestKillBFunction(InputUnsafe inp, Config config, float intthresh)
		{
			_TestKillBFunction(inp.dataInt, inp.barr, intthresh, config.lengthmz, config.numz, config.isolength, inp.isotopeops, inp.isotopeval); 
		}
	}
}
