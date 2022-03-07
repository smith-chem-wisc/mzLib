using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using System.Runtime.InteropServices; 

namespace TestUniDec
{
	class TestHomework
	{
		// This class shows how to use methods from unmanaged DLLs, i.e. an unmanaged class API.  
		// The template is 
		/*	
		 *	[DllImport()]
		 *	private static extern type _DllFunction(); 
		 *	public static type DllFunction(); 
		 *	
		 *	A public method should be paired with the private method because in some 
		 *	instances, we will want to do error handling and or memory allocation/deallocation in the public method 
		 *	in order to have neater API code. 
		 *	
		 * General guidelines: 
		 * 
		 * Build the C .dll in Release mode. You need to add __declspec(dllexport) to all function you are planning to export. 
		 * 
		 * Copy and paste the C generated .dll into the same repository as your API class, then go to the .dll Properties and change to Copy Always. 
		 * 
		 * Turn off pre-compiled headers in the C code. It won't compile with them.
		 * 
		 * Change all files to .c extension. Otherwise, you'll compile C++, and that requires different skill sets. 
		 * 
		 * Don't bother trying to figure out how to use C++ calls in your unmanaged code. It's not worth it. Just write C code. 
		 * 
		 * Try not to use unsafe methods. If you do want to use unsafe methods, wrap it in an unsafe block. Read up on the fixed keyword 
		 * for pointer allocation in C#. 
		 * 
		 * In general, you want to perform very little memory allocation in your unmanaged native code. 
		 * 
		 * If you perform memory allocation in C, you must also pass a method to the API that frees the memory you passed. 
		 * Name this method with a _dispose suffix. Call it in the public method. Example: 
		 * [DllImport()] 
		 * private static extern ReturnStruct _MethodThatAllocatesMemory(parameters); 
		 * 
		 * [DllImport()] 
		 * private static extern void _MethodThatRemoveObjectFromMemory(*pointerToAllocatedMemory); 
		 * 
		 * public static void MethodThatAllocatesMemory(parameters, out ReturnStruct result){
		 *		// Note that if a pointer is a component of ReturnStruct, you'll need to process the pointer
		 *		// before removing the object from memory. Otherwise, you'll return nothing and probably throw an error. 
		 *		result = _MethodThatAllocatesMemory(parameters); 
		 *		_MethodThatRemovesObjectFromMemory(*pointerToAllocatedMemory);			
		 * }
		 * 
		 * Good luck. 
		 */

		[DllImport("UniDecFunctions", EntryPoint = "PrintHelloWorld")]
		private static extern void printHelloWorld(); 
		public static void PrintHelloWorld()
		{
			printHelloWorld(); 
		}
		[DllImport("UniDecFunctions", EntryPoint = "ReturnInteger")]
		private static extern int returnInteger(); 
		public static int ReturnInteger()
		{
			return returnInteger(); 
		}

		[DllImport("UniDecFunctions", EntryPoint = "AddTwoIntegers")]
		private static extern int _AddTwoIntegers(int i, int j);
		public static int AddTwoIntegers(int i, int j)
		{
			// need to validate that the unmanaged code returned the expected result. 
			// ...eventually... 
			int result = _AddTwoIntegers(i, j);
			return result; 
		}

		[DllImport("UniDecFunctions", EntryPoint = "MultiplyTwoDouble")]
		private static extern double _MultiplyTwoDouble(double i, double j);
		public static double MultiplyTwoDouble(double i, double j)
		{
			double result = _MultiplyTwoDouble(i, j);
			return result; 
		}
		[DllImport("UniDecFunctions", EntryPoint = "EchoMessage")]
		// Pointers and strings passing: https://www.codeproject.com/Articles/1189085/Passing-strings-between-managed-and-unmanaged-code
		private static extern IntPtr _EchoMessage([MarshalAs(UnmanagedType.LPWStr )] string message); 
		public static string EchoMessage(string message)
		{
			IntPtr ptr = _EchoMessage(message);
			string echoedMessage = Marshal.PtrToStringUni(ptr);
			return echoedMessage; 
			
		}

		[DllImport("UniDecFunctions.dll", EntryPoint = "PassIntArrayModifyThenReturnArray")]
		private static extern IntPtr _PassIntArrayModifyThenReturnArray(int[] arrayOfInt, int length); 
		public static int[] PassIntArrayModifyThenReturnArray(int[] arrayOfInt, int length)
		{
			IntPtr ptrResults = _PassIntArrayModifyThenReturnArray(arrayOfInt, length);
			int[] result = new int[length];  
			Marshal.Copy(ptrResults, result, 0, length);
			return result; 
		}

		// Next functions show how marshaling (should) work. 
		// See https://docs.microsoft.com/en-us/dotnet/standard/native-interop/type-marshaling
		// For C code: https://docs.microsoft.com/en-us/dotnet/framework/interop/marshaling-data-with-platform-invoke
		// for a list of what C and C++ types correspond to which C# types. 
		[DllImport("UniDecFunctions.dll", EntryPoint = "InitializeAndReturnStruct")]
		private static extern TestStruct _InitializeAndReturnStruct(); 
		public static TestStruct InitializeAndReturnStruct()
		{
			TestStruct results = _InitializeAndReturnStruct();
			return results; 
		}
		[DllImport("UniDecFunctions.dll", EntryPoint = "SendingAMessage")]
		private static extern string _SendingAMessage(); 
		public static string SendingAMessage()
		{
			string result = _SendingAMessage();
			return result;
		}
		[DllImport("UniDecFunctions.dll", EntryPoint = "InitTestStruct2")]
		private static extern TestStruct2 _InitTestStruct2(int val1, int[] array1);

		public static TestStruct2 InitTestStruct2(int val1, int[] array1)
		{
			// if you are passing a pointer from C back to C#, the corresponding struct that 
			// the C code returns to needs contain an IntPtr, not the final array.
			// Then you need to use Marshal.Copy() to conver the IntPtr to the correct array. 
			TestStruct2 ts = _InitTestStruct2(val1, array1);
			return ts; 
		}
		[DllImport("UniDecFunctions.dll", EntryPoint = "InitTestStruct2")]
		private static extern TestStruct2Unsafe _InitTestStruct2ToUnsafe(int val1, int[] array1);
		public static TestStruct2Unsafe InitTestStruct2ToUnsafe(int val1, int[] array1)
		{
			TestStruct2Unsafe ts2Unsafe = _InitTestStruct2ToUnsafe(val1, array1);
			return ts2Unsafe; 
		}
	}
	class TestHomeworkDLLFunctions
	{
		[Test]
		public void TestPrintTestFunction()
		{
			TestHomework.PrintHelloWorld();
		}
		[Test]
		public void TestReturnInteger()
		{
			int result = TestHomework.ReturnInteger();
			Assert.AreEqual(1, result);
		}
		[Test]
		public void TestPassIntArrayModifyThenReturnArray()
		{
			int[] testArray = { 0, 1, 2, 3, 4, 5 };
			int[] results = TestHomework.PassIntArrayModifyThenReturnArray(testArray, testArray.Length);
			foreach (int i in results)
			{
				TestContext.WriteLine(i.ToString());
			}
		}
		[Test]
		public void TestEchoMessage()
		{
			string s = "This is echoed from C.";
			string results = TestHomework.EchoMessage(s);
			TestContext.WriteLine(results);
		}
		[Test]
		public void TestTestStruct()
		{
			TestStruct tsRes = TestHomework.InitializeAndReturnStruct();
			Console.WriteLine(tsRes.val1.ToString());
		}
		[Test]
		public void TestSendingAMessage()
		{
			string result = TestHomework.SendingAMessage();
			Console.WriteLine(result);
		}
		[Test]
		public void TestTestStruct2()
		{
			int val1 = 27;
			int[] array1 = { 1, 2, 3, 4, 5 };
			TestStruct2 ts = TestHomework.InitTestStruct2(val1, array1);
			Console.WriteLine(ts.val1.ToString());
			Console.WriteLine(string.Join("; ", ts.array1));
		}
		[Test]
		public void TestTestStruct2ReturnIntArray()
		{
			int val1 = 2;
			int[] array1 = { 1, 2, 3, 4, 5 };
			// Idiom that's important to know. 
			int[] res = new int[array1.Length];
			TestStruct2 ts = TestHomework.InitTestStruct2(val1, array1);
			Marshal.Copy(ts.array1, res, 0, res.Length);
			Console.WriteLine(string.Join("; ", res));
		}
		[Test]
		public void TestConvertStruct2ToClass()
		{
			/*
			 * Shows how to convert a struct that matches native, unmanaged code to a class. 
			 * 1) Return the initial struct from the native code. 
			 * 2) Create a matching class that takes the struct as an input argument in its constructor. 
			 */
			int val1 = 2;
			int[] array1 = { 1, 2, 3, 4, 5 };
			int[] res = new int[array1.Length];
			TestStruct2 ts = TestHomework.InitTestStruct2(val1, array1);
			Struct2 struct2 = new Struct2(ts, array1.Length);
			Assert.AreEqual(array1.Length, struct2.Array1.Length);
			Assert.AreEqual(array1, struct2.Array1); 
		}
		[Test]
		public void AccessTestStruct2Unsafe()
		{
			// This method shows how you can use unsafe code to directly return pointers from 
			// unmanaged code as well. I'm not sure what would be easier to do... 
			// I think IntPtr would be more safe, but as long as evreything is done on the stack, 
			// I don't see why it would hurt to use unsafe code in C#. 
			int val1 = 20;
			int[] array1 = { 1, 2, 3, 4, 5 };
			unsafe
			{
				TestStruct2Unsafe ts2Unsafe = TestHomework.InitTestStruct2ToUnsafe(val1, array1);
				for(int i = 0; i < array1.Length; i++)
				{
					Console.WriteLine(ts2Unsafe.array1[i]); 
				}
				// You need to be careful about the memory assignment. 
				/*	
				 * Note the types we are using. In this case we have an array of int[]. 
				 * This type is allocated on the stack, not the heap. Memory allocated to the stack 
				 * gets erased once we lose the local scope. In this case, we can be sure that we don't 
				 * need to explicity get rid of the memory allocatd to array1. 
				 * 
				 * Structs are also allocated on the stack, so that means that they should also be deleted once 
				 * removed from the function scope. 
				 * https://stackoverflow.com/questions/4853213/are-structs-always-stack-allocated-or-sometimes-heap-allocated#:~:text=your%20specific%20question-,Are%20struct%20instances%20sometimes%20allocated%20on%20the%20heap%3F,sometimes%20allocated%20on%20the%20heap.
				 */
				fixed (int* ptr = &array1[0])
				{
					// pointers allocated in fixed statements will be handled by the garbage collector at the end 
					// of the fixed statement. 

					// Console prints out that we are using the same addres in ts2Unsafe.array1 as we were when 
					// we initialized array1. Pretty neat. 
					Console.WriteLine(string.Join("; ",(IntPtr)ptr, (IntPtr)ts2Unsafe.array1)); 

					
				}
				// Theoretically at the end of the function, the garbage collector should clean up array1, thus keeping 
				// us memory safe. 
			}
		}
	}
}
