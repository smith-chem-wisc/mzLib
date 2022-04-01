using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class MassIntensityDetermination
	{
		public static void IntegrateMassIntensities(Config config, Decon decon, InputUnsafe inp)
		{
			float massmax = config.masslb;
			float massmin = config.massub;
			fixed (float* massaxisPtr = &decon.massaxis[0],
				massaxisvalPtr = &decon.massaxisval[0],
				blurPtr = &decon.blur[0],
				massgridPtr = &decon.massgrid[0],
				newblurPtr = &decon.newblur[0])
			{
				if (config.poolflag == 0)
				{
					if (config.rawflag == 1 || config.rawflag == 3)
					{
						IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin,
							decon.mlen, massaxisPtr, massaxisvalPtr, blurPtr, massgridPtr);
					}
					if (config.rawflag == 0 || config.rawflag == 2)
					{
						IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen,
							 massaxisPtr, massaxisvalPtr, newblurPtr, massgridPtr);
					}
				}
				else if (config.poolflag == 1)
				{
					if (config.rawflag == 1 || config.rawflag == 3)
					{
						InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, massaxisPtr,
							config.adductmass, inp.dataMZ, massgridPtr, massaxisvalPtr, blurPtr);
					}
					if (config.rawflag == 0 || config.rawflag == 2)
					{
						InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab,
							massaxisPtr, config.adductmass, inp.dataMZ, massgridPtr, massaxisvalPtr, newblurPtr);
					}
				}
				else if (config.poolflag == 2)
				{
					if (config.rawflag == 1 || config.rawflag == 3)
					{
						SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, massaxisPtr, config.adductmass,
							inp.dataMZ, massgridPtr, massaxisvalPtr, blurPtr);
					}
					if (config.rawflag == 0 || config.rawflag == 2)
					{
						SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, massaxisPtr, config.adductmass,
							inp.dataMZ, massgridPtr, massaxisvalPtr, newblurPtr);
					}
				}
				else
				{
					// throw new error
				}
			}
		}
		[DllImport("TestDLL.dll", EntryPoint = "IntegrateTransform")]
		private static extern void _IntegrateTransform(int lengthmz, int numz, float* mtab, float massmax, float massmin, int maaxle,
			float* massaxis, float* massaxisval, float* blur, float* massgrid);
		public static void IntegrateTransform(int lengthmz, int numz, float* mtab, float massmax,
			float massmin, int maaxle, float* massaxis, float* massaxisval, float* blur, float* massgrid)
		{
			_IntegrateTransform(lengthmz, numz, mtab, massmax, massmin, maaxle,
			massaxis, massaxisval, blur, massgrid);
		}
		[DllImport("TestDLL.dll", EntryPoint = "InterpolateTransform")]
		private static extern void _InterpolateTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
			float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur);
		public static void InterpolateTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
			float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur)
		{
			_InterpolateTransform(maaxle, numz, lengthmz, nztab, massaxis, adductmass,
				dataMZ, massgrid, massaxisval, blur);
		}
		[DllImport("TestDLL.dll", EntryPoint = "SmartTransform")]
		private static extern void _SmartTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
			float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur);
		public static void SmartTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
			float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur)
		{
			_SmartTransform(maaxle, numz, lengthmz, nztab, massaxis, adductmass,
				dataMZ, massgrid, massaxisval, blur);
		}
	}
}
