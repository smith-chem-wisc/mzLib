﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public static class APITesting
	{
		[DllImport("TestDLL.dll", EntryPoint = "Average")]
		private static extern float Average(int length, float[] xarray);
		public static float Average(float[] xarray) 
		{
			return Average(xarray.Length, xarray);
		}
		[DllImport("TestDLL.dll", EntryPoint = "ExecuteTestFFT")]
		private static extern int _ExecuteTestFFT();
		public static int ExecuteTestFFT()
		{
			return _ExecuteTestFFT(); 
		}
		[DllImport("TestDLL.dll", EntryPoint = "max1d")]
		private static extern float _max1d(float[] blur, int lengthmz); 
		public static float Max1d(float[] blur)
		{
			return _max1d(blur, blur.Length); 
		}
	}
}
