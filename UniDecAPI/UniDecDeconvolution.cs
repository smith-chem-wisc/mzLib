using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static class UniDecDeconvolution
	{ 
		public static void DeconvoluteWithUniDec(this MzSpectrum spectrum, out DeconResults deconResults)
		{
			UniDecDeconEngine deconEngine = new();
			double[] yarray = deconEngine.LinearizeMassSpectrum(spectrum, 0.5, out double[] xarray);
			MzSpectrum linSpectrum = new(xarray, yarray, false); 
			deconResults = deconEngine.Run(spectrum); 
		}
	}
	public class UniDecDeconEngine 
	{
		Config config;
		DeconUnsafe decon;
		InputUnsafe inp;
		UnmanagedHandling unHandler;
		byte[] barr; 

		public DeconResults Run(MzSpectrum spectrum)
		{
			DeconResults results; 
			using (unHandler = new UnmanagedHandling())
			{
				Setup(spectrum);
				Processing();
				results = ScoreAndOutput(); 
			}
			return results; 
		}
		private void RawDataPreprocessing(MzSpectrum spectrum)
		{
			double minmz, maxmz, gaussianSmoothWidth, binSize, intthresh, dataReductionPercent;
			bool baseLineSubtration, smooth; 

			// crop data

			// Smooth Data
			// Linearize Data
			// BaseLineSubtraction
			// DataReduction
			// Normalization
			// Intensity Thresholding
			// DataReduction
			// ZeroRemoval
			// Duplicate Removal. 
		}
		// Trims m/z range of data. 
		public MzSpectrum DataChop(MzSpectrum spectrum, double minMz, double maxMz)
		{
			return new MzSpectrum(spectrum.Extract(minMz, maxMz)); 
		}
		public double[] GaussianSmoothingFilter(MzSpectrum spectrum, double sigma, int kernelSize)
		{
			double deviation = (kernelSize - 1) / 6;
			double[] kernel = GenerateGuassianKernel(deviation, kernelSize);
			return ConvoluteFull(spectrum.YArray, kernel); 
		}
		public double[] ConvoluteFull(double[] signal, double[] kernel)
		{
			// Output needs padding. Length shouuld be 
			double[] output = new double[signal.Length + 2 * kernel.Length - 1];
			// create zero padded array 
			Array.Copy(signal, 0, output, kernel.Length - 1, signal.Length); 
			// outer loop iterates through the signal
			double[] accumulator = new double[output.Length]; 
			
			for(int sig = kernel.Length; sig < output.Length - kernel.Length; sig++)
			{
				double sum = 0; 
				// inner loop iterates through the kernel
				// kernel values need to go from 
				for(int kern = 0; kern < kernel.Length; kern++)
				{
					accumulator[sig] += kernel[kern] * output[sig - kern]; 
				}
			}

			return RemovePadding(accumulator, signal.Length); 
		}
		private double[] RemovePadding(double[] array, int originalDimensions)
		{
			int padding = array.Length - originalDimensions; 
			List<double> arrayList = array.ToList();
			return arrayList.Skip(padding / 2).Take(originalDimensions).ToArray(); 
		}
		public void ConvoluteFull(ref double[] signal, double[] kernel)
		{
			signal = ConvoluteFull(signal, kernel); 
		}
		public double[,] ConvoluteFull2D(double[,] matrix, double[,] kernel)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)]; 
			for(int y = 0; y < matrix.GetLength(1); y++)
			{
				for(int x = 0; x < matrix.GetLength(0); x++){
					double sum = 0; 

					for(int kernelY = -kernel.GetLength(0)/2; kernelY < kernel.GetLength(0)/2; kernelY++)
					{
						for(int kernelX = -kernel.GetLength(1)/2; kernelX < kernel.GetLength(1)/2; kernelX++)
						{
							int sourceY = y + kernelY;
							int sourceX = x + kernelX;

							if (sourceX < 0)
								sourceX = 0;

							if (sourceX >= matrix.GetLength(0))
								sourceX = matrix.GetLength(0) - 1;

							if (sourceY < 0)
								sourceY = 0;

							if (sourceY >= matrix.GetLength(1))
								sourceY = matrix.GetLength(1) - 1;

							sum += matrix[sourceX, sourceY]; 
						}
					}
					result[x, y] = sum; 
				}
			}
			
			return result;
		}
		public double[] GenerateGuassianKernel(double sigma, int size)
		{
			/* Gaussian formula: 
			 * g(x) = alpha * exp{-(x-mu)^2/(2 * sigma^2)}
			 * where alpha = 1/(sigma * sqrt(2*pi)) 
			 */

			double[] kernel = new double[size];
			int half = size / 2;

			double alpha = 1 / (sigma * Math.Sqrt(2 * Math.PI)); 

			for(int i = 0; i < size; i++)
			{
				double beta = -(i - half) * (i - half) / (2 * sigma * sigma);
				kernel[i] = alpha * Math.Exp(beta);
			}
			return kernel; 
		}
		public double[,] GenerateGaussianKernel2D(double sigma, int size)
		{
			double[,] kernel = new double[size, size];
			double frontTerm = 1 / (2 * Math.PI * sigma * sigma);
			double twoSigSquared = 2 * sigma * sigma; 

			for(int i = 0; i < kernel.GetLength(0); i++)
			{
				for(int j = 0; j < kernel.GetLength(1); j++)
				{
					kernel[i, j] = frontTerm * Math.Exp(-(i * i + j * j) / twoSigSquared); 
				}
			}
			return kernel; 

		}
		private void LinearizeData() { }
		private void NonLinearizeData() { }
		private void RemoveNoise() { }
		private void Normalize()
		{

		}
		private void IntensityThresholding() { }
		
		private void PeakDetect()
		{

		}
		private void RemoveMiddleZeroes()
		{

		}
		private void RemoveDuplicateValues() { }
		private void Setup(MzSpectrum spectrum)
		{
			unsafe
			{
				// create the config, input and decon structs. 
				UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(spectrum, out config);
				inp = UniDecAPIMethods.InputMethods.SetupInputs();
				decon = new DeconUnsafe();

				// copy the x- and y-axis to a pointer and assign it to the input struct. 
				inp.dataInt = (float*)unHandler.AllocateToUnmanaged(spectrum.YArray.ConvertDoubleArrayToFloat());
				inp.dataMZ = (float*)unHandler.AllocateToUnmanaged(spectrum.XArray.ConvertDoubleArrayToFloat());
				
				// further inp setup after the data from the spectrum is filled. 
				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
				UniDecAPIMethods.UtilityMethods.SetLimits(config, inp);
				IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
				config.isolength = Math.Abs(isoStruct.isolength);

				inp.isotopeops = (int*)unHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(int));
				inp.isotopeval = (float*)unHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(float));

				Isotopes.MakeIsotopesFromAPI(config, inp, isoStruct);

				// I believe barr stands for "Bad array". 0 if bad; 1 if good. 
				int numberOfElementsInBarr = config.lengthmz * config.numz;

				barr = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberOfElementsInBarr);
				// C is returning a char, so this function needs to convert the char to the equivalent C# byte. 
				UniDecAPIMethods.UtilityMethods.ConvertASCIIBytesFromCToByteInCS(ref barr);
			}
		}
		private unsafe void Processing()
		{
			float threshold = config.psthresh * Math.Abs(config.mzsig) * config.peakshapeinflate;

			int* starttab = (int*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));
			int* endtab = (int*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));

			int maxlength = Convolution.SetStartEnds(config, ref inp, starttab, endtab, threshold);

			int pslen = config.lengthmz * maxlength;
			float* mzdist = (float*)unHandler.AllocateToUnmanaged(pslen, typeof(float));
			float* rmzdist = (float*)unHandler.AllocateToUnmanaged(pslen, typeof(int));
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int* mind = (int*)unHandler.AllocateToUnmanaged(mlength, typeof(int));
			float* mdist = (float*)unHandler.AllocateToUnmanaged(mlength, typeof(int));
			for (int i = 0; i < mlength; i++)
			{
				mind[i] = i - (mlength - 1) / 2;
				if (config.msig != 0)
				{
					mdist[i] = (float)(Math.Exp(-Math.Pow((i - (zlength - 1) / 2.0), 2)) / (2.0 * config.zsig * config.zsig));
				}
				else
				{
					mdist[i] = 1;
				}
			}
			int* zind = stackalloc int[zlength];
			float* zdist = stackalloc float[zlength];

			int numclose = mlength * zlength;
			int* closemind = (int*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closezind = (int*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			float* closeval = (float*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closeind = (int*)unHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));
			float* closearray = (float*)unHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));

			for (int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
			}

			Normalization.simp_norm_sum(mlength, mdist);
			Normalization.simp_norm_sum(zlength, zdist);
			Normalization.simp_norm_sum(numclose, closeval);

			Blur._MakeSparseBlur(numclose, inp.barr, closezind, closemind, inp.mtab, inp.nztab,
				inp.dataMZ, closeind, closeval, closearray, config);

			int badness = 1;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					badness = 0;
				}
			}
			if (badness == 1)
			{
				throw new InvalidOperationException("Badness = 1...");
			}

			float dmax = MathUtilities.Max(inp.dataInt, config.lengthmz);
			float betafactor = 1;
			if (dmax > 1)
			{
				betafactor = dmax;
			}

			UniDecAPIMethods.UtilityMethods.KillB_CS(inp.dataInt, ref barr, config.intthresh, config.lengthmz, config.numz, config.isolength,
				inp.isotopeops, inp.isotopeval);
			// DirectUniDecPort.FitFunctions.KillB(inp.dataInt, barr, 
			//	config.intthresh, config.lengthmz, config.numz, 
			//	config.isolength, inp.isotopeops, inp.isotopeval); 

			decon.blur = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			decon.newblur = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			float* oldblur = stackalloc float[config.lengthmz * config.numz];
			decon.baseline = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			decon.noise = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));

			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / ((float)(config.numz + 2));
				if (config.baselineflag == 1)
				{

					decon.baseline[i] = val;
					decon.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						if (config.isotopemode == 0)
						{
							decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = val;
						}
						else { decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
					}
					else
					{
						decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 0;
					}
				}
			}

			// copy decon.blur to oldblur: 
			oldblur = decon.blur;
			// copy decon.newblur to decon.blur: 
			decon.newblur = decon.blur;

			float[] dataInt2 = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.dataInt, config.lengthmz);


			Convolution.DeconvolveBaseline(config, inp, decon);
			float blurmax = 0F;
			decon.conv = 0F;

			fixed (float* dataInt2Ptr = &dataInt2[0])
			{
				Blur.PerformIterations(ref decon, config, inp, betafactor, maxlength,
					starttab, endtab, mzdist, numclose, closeind, closearray, zlength, mlength,
					closemind, closezind, mdist, dataInt2Ptr, zdist, barr, rmzdist, oldblur);
			}


			if (config.peakshapeinflate != 1 && config.mzsig != 0)
			{
				if (config.speedyflag == 0)
				{
					MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse: 0, starttab, endtab, maxlength);
				}
				else
				{
					MZPeak.MakePeakShape1D(config, inp, threshold, mzdist, rmzdist, makeReverse: 0);
				}
			}

			blurmax = MathUtilities.Max(decon.blur, config.lengthmz * config.numz);
			float cutoff = 0F;
			if (blurmax != 0)
			{
				cutoff = 0.000001F;
			}

			ArrayIndexing.ApplyCutoff1D(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);

			decon.fitdat = (float*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(float));


			decon.error = ErrorFunctions.ErrFunctionSpeedy(config, ref decon, barr, inp.dataInt, maxlength, inp.isotopeops, inp.isotopeval,
				starttab, endtab, mzdist);
			//decon.error = ErrorFunctions.ErrFunctionSpeedy(config, decon, barr, inp.dataInt, maxlength,
			//	inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist);

			if (config.intthresh != -1)
			{
				for (int i = 0; i < config.lengthmz - 1; i++)
				{
					if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
					{
						decon.fitdat[i] = 0F;
						decon.fitdat[i + 1] = 0F;
					}
				}
			}
			// not tested yet. 
			if (config.isotopemode == 2)
			{
				Isotopes.MonoisotopicToAverageMass(config, inp, decon, barr);
			}

			float newblurmax = blurmax;
			if (config.rawflag == 0 || config.rawflag == 2)
			{
				if (config.mzsig != 0)
				{
					newblurmax = Convolution._Reconvolve(config.lengthmz, config.numz, maxlength,
						starttab, endtab, mzdist, decon.blur, decon.newblur, config.speedyflag, inp.barr);
				}
				else
				{
					decon.newblur = decon.blur;
				}
			}
			float massmax = config.masslb;
			float massmin = config.massub;

			if (config.fixedmassaxis == 0)
			{
				for (int i = 0; i < config.lengthmz; i++)
				{
					for (int j = 0; j < config.numz; j++)
					{
						if (decon.newblur[ArrayIndexing.Index2D(config.numz, i, j)] * barr[ArrayIndexing.Index2D(config.numz, i, j)] > newblurmax * cutoff)
						{
							float testmax = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] + threshold * inp.nztab[j] + config.massbins;
							float testmin = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] - threshold * inp.nztab[j];

							//To prevent really wierd decimals
							testmin = (float)Math.Round(testmin / config.massbins) * config.massbins;
							testmax = (float)Math.Round(testmax / config.massbins) * config.massbins;

							if (testmax > massmax) { massmax = testmax; }
							if (testmin < massmin) { massmin = testmin; }
						}
					}
				}
			}
			else { massmax = config.massub; massmin = config.masslb; }

			//Checks to make sure the mass axis is good and makes a dummy axis if not
			decon.mlen = (int)((massmax - massmin) / config.massbins);
			if (decon.mlen < 1)
			{
				massmax = config.massub;
				massmin = config.masslb;
				decon.mlen = (int)((massmax - massmin) / config.massbins);

				//Declare the memory

				decon.massaxis = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massaxisval = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massgrid = (float*)unHandler.AllocateToUnmanaged(decon.mlen * config.numz, typeof(float));

				//Create the mass axis
				for (int i = 0; i < decon.mlen; i++)
				{
					decon.massaxis[i] = massmin + i * config.massbins;
				}
				decon.uniscore = 0;
			}
			else
			{

				//Declare the memory
				decon.massaxis = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massaxisval = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massgrid = (float*)unHandler.AllocateToUnmanaged(decon.mlen * config.numz, typeof(float));


				//Create the mass axis
				for (int i = 0; i < decon.mlen; i++)
				{
					decon.massaxis[i] = massmin + i * config.massbins;
				}
			}

			MassIntensityDetermination.IntegrateMassIntensities(config, ref decon, inp);
		}
		private unsafe DeconResults ScoreAndOutput()
		{
			float scorethresh = 0f;
			config.peakwin = 20;

			decon.peakx = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
			decon.peaky = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));

			decon.uniscore = Scoring.UniScorePorted(config, ref decon, inp, scorethresh, config.peakwin, unHandler);

			return new DeconResults(decon, config); 
		}
		public double[,] CreateChargeAndMZMatrix(double[] intArray, int minCharge, int maxCharge)
		{
			// create charge range
			int[] chargeArray = Enumerable.Range(minCharge, maxCharge - minCharge).ToArray();
			double[,] chargeAndMzMatrix = new double[chargeArray.Length, intArray.Length]; 

			for(int i = 0; i < chargeAndMzMatrix.GetLength(0); i++)
			{
				for(int j = 0; j < chargeAndMzMatrix.GetLength(1); j++)
				{
					chargeAndMzMatrix[i, j] = chargeArray[i] * intArray[j];
				}
			}
			return chargeAndMzMatrix; 
		}
		public double[] SumChargeAxis(double[,] chargeMzMatrix)
		{
			double[] colSums = new double[chargeMzMatrix.GetLength(1)]; 
			for(int i = 0; i < colSums.Length; i++)
			{
				colSums[i] = Enumerable.Range(0, chargeMzMatrix.GetLength(0)).Select(x => chargeMzMatrix[x, i]).Sum();
			}
			return colSums; 
		}
		public double[] ConvoluteSummedChargeAxis(double[] chargeAxis, double[] peakShapeKernel)
		{
			return ConvoluteFull(chargeAxis, peakShapeKernel); 
		}
		public double[,] SmoothDeltaFunctions(double[,] matrix, double[,] kernel)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
			//convolute the kernel with the 2D image. 
			return ConvoluteFull2D(matrix, kernel); 

		}
		public void SmoothDeltaFunctions(ref double[,] matrix, double[,] kernel)
		{
			matrix = SmoothDeltaFunctions(matrix, kernel); 
		}
		public double[,] MatrixMultiplication(double[,] matrix1, double[,] matrix2)
		{
			int rowsMatrix1 = matrix1.GetLength(0);
			int colsMatrix1 = matrix1.GetLength(1);
			int rowsMatrix2 = matrix2.GetLength(0);
			int colsMatrix2 = matrix2.GetLength(1); 
			double[,] result = new double[rowsMatrix1, colsMatrix2];
			for(int i = 0; i < colsMatrix2; i++)
			{
				for(int j = 0; j < rowsMatrix1; j++)
				{
					double sum = 0; 
					for(int k = 0; k < colsMatrix1; k++)
					{
						sum += matrix1[j, k] * matrix2[k, i]; 
					}
					result[j, i] = sum; 
				}
			}
			return result; 
		}

		public double[] MultiplyArrayByScalar(double[] array, double scalar)
		{
			return array.Select(i => i * scalar).ToArray(); 
		}
		public double[] TakeReciprocalOfArray(double[] array)
		{
			return array.Select(i => i = (i != 0) ? 1 / i : 0).ToArray(); 
		}
		public double[,] ElementwiseMultiplication(double[,] matrix1, double[,] matrix2)
		{
			double[,] result = new double[matrix1.GetLength(0), matrix1.GetLength(1)]; 

			for(int i = 0; i < matrix1.GetLength(0); i++)
			{
				for(int j = 0; j < matrix1.GetLength(1); j++)
				{
					result[i, j] = matrix1[i,j] * matrix2[i,j]; 
				}
			}
			return result; 
		}
		public double[,] ElementwiseMutiplyArray(double[,] matrix, double[] array)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)]; 
			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				for(int j = 0; j < matrix.GetLength(1); j++)
				{
					result[i, j] = matrix[i, j] * array[j]; 
				}
			}
			return result; 
		}
		public double[] ElementwiseMultiplyArrays(double[] array1, double[] array2)
		{
			double[] result = new double[array1.Length]; 
			for(int i = 0; i < result.Length; i++)
			{
				result[i] = array1[i] * array2[i]; 
			}
			return result; 
		}
		public double[,] CreateDeepCopy(double[,] source) 
		{
			double[,] copy = new double[source.GetLength(0), source.GetLength(1)];
			double[] buffer = new double[source.GetLength(0) * source.GetLength(1)];

			Buffer.BlockCopy(source, 0, buffer, 0, buffer.Length * sizeof(double));
			Buffer.BlockCopy(buffer, 0, copy, 0, buffer.Length * sizeof(double));
			return copy; 
		}
		public double[,] CreateChargeByMzTable(int[] charges, double[] mzAxis, double adductMass)
		{
			double[,] massTable = new double[charges.Length, mzAxis.Length]; 
			for(int i = 0; i < massTable.GetLength(0); i++)
			{
				for(int j = 0; j < massTable.GetLength(1); j++)
				{
					massTable[i, j] = charges[i] * mzAxis[j] - adductMass * charges[i];
				}
			}
			return massTable; 
		}
		public double[,] CreateDeconvolutedMassTable(double[] massAxis, int[] chargeArray)
		{
			return new double[chargeArray.Length, massAxis.Length];
		}
		public void IntegrateTransform(double[,] massTable, double[] massaxis, ref double[] massaxisVal, 
			double[,] ft, ref double[,] deconMassTable)
		{
			// this function uses an accumulator grid of charge vs mass. 
			// The accumulator is the deconMassTable. It is made up of the massaxis on the x-axis
			// (cols) and the charge on the y-axis (rows). 

			// if a mass from the mass table matches a mass from the deconvoluted mass table, 
			// the value of that is pulled from the ft grid and added to the accumulator grid. 
			// the location of the addition is the (charge, index) where index = the matching index 
			// between the mass from the mass index and the mass axis. 
			double massMax = massaxis.Max();
			double massMin = massaxis.Min(); 
			// i is charge, j is mass
			for(int i = 0; i < massTable.GetLength(0); i++)
			{
				for(int j = 0; j < massTable.GetLength(1); j++)
				{
					double testmass = massTable[i,j];
					if (testmass > massMax || testmass < massMin) continue; 

					int index = Array.BinarySearch(massaxis, testmass);
					if (index < 0) index = ~index;
					
					double newval = ft[i, j];
					if(massaxis[index] == testmass)
					{
						massaxisVal[index] += newval;
						deconMassTable[i,index] += newval; 
					}
					if(massaxis[index] < testmass && index < massaxis.Length - 2)
					{
						int index2 = index + 1;
						double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
						massaxisVal[index] += (1.0 - interpos) * newval;
						deconMassTable[i,index] += (1.0 - interpos) * newval;

						massaxisVal[index2] += interpos * newval;
						deconMassTable[i,index2] += interpos * newval; 
					}
					if(index > 0 && massaxis[index] > testmass)
					{
						int index2 = index - 1;
						double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
						massaxisVal[index] += (1 - interpos) * newval;
						deconMassTable[i,index] += (1.0 - interpos) * newval;
						massaxisVal[index2] += interpos * newval;
						deconMassTable[i,index2] += interpos * newval;
					}
				}
			}
			
		}
		public double LinearInterpolatePosition(double x1, double x2, double x)
		{
			if(x2 - x1 == 0)
			{
				return 0; 
			}
			return (x - x1) / (x2 - x1); 
		}
		public void CreateMassAxes(out double[] massAxis, out double[] massAxisVals, 
			double massMax, double massMin, double massBinWidth)
		{
			double massDiff = massMax - massMin;
			int elementsOfMassAxis = (int)(massDiff / massBinWidth);

			massAxis = new double[elementsOfMassAxis];
			massAxisVals = new double[elementsOfMassAxis]; 

			for(int i = 0; i < elementsOfMassAxis; i++)
			{
				massAxis[i] = massMin + i * massBinWidth; 
			}
		}
		public double[] RowSums(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(0)];
			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				result[i] = Enumerable.Range(0, matrix.GetLength(1))
					.Select(x => matrix[i,x])
					.Sum();
			}
			return result;
		}
		public double[] ColSums(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(1)]; 
			for(int i = 0; i < matrix.GetLength(1); i++)
			{
				result[i] = Enumerable.Range(0, matrix.GetLength(0))
					.Select(x => matrix[x, i])
					.Sum(); 
			}
			return result; 
		}
		public double[,] ApplyLogMeanFilter(double[,] matrix, int width)
		{
			double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
			double frontTerm = 1 / (2 * width + 1); 

			for (int y = 0; y < matrix.GetLength(1); y++)
			{
				for (int x = 0; x < matrix.GetLength(0); x++)
				{
					double sum = 0;
					for (int i = -width/2; i < width/2; i++)
					{
						int sourceY = y + i;
						int sourceX = x + i;

						if (sourceX < 0)
							sourceX = 0;

						if (sourceX >= matrix.GetLength(0))
							sourceX = matrix.GetLength(0) - 1;

						if (sourceY < 0)
							sourceY = 0;

						if (sourceY >= matrix.GetLength(1))
							sourceY = matrix.GetLength(1) - 1;

						if(matrix[sourceX, sourceY] > 0)
						{
							sum += Math.Log(matrix[sourceX, sourceY]);
						}
					}
					result[x, y] = Math.Exp(sum * frontTerm);
				}
			}
			return result;
		}
		public double[] ApplyLogMeanFilter(double[] array, int width)
		{
			double[] result = new double[array.Length];
			double frontTerm = 1 / (2 * width + 1);

			
				for (int x = 0; x < array.Length; x++)
				{
					double sum = 0;
					for (int i = -width / 2; i < width / 2; i++)
					{
						int sourceX = x + i;

						if (sourceX < 0)
							sourceX = 0;

						if (sourceX >= array.GetLength(0))
							sourceX = array.GetLength(0) - 1;


						if (array[sourceX] > 0)
						{
							sum += Math.Log(array[sourceX]);
						}
					}
					result[x] = Math.Exp(sum * frontTerm);
				}
			return result;
		}
		public double[] LinearizeMassSpectrum(MzSpectrum spectrum, double binWidth, out double[] mzAxisNew)
		{
			double firstpoint = Math.Ceiling(spectrum.XArray[0] / binWidth) * binWidth;
			double lastpoint = Math.Floor(spectrum.XArray.Last() / binWidth) * binWidth;
			int numberBins = (int)((lastpoint - firstpoint) / binWidth);
			mzAxisNew = Enumerable.Range(0, numberBins).Select(i => firstpoint + i * binWidth).ToArray();

			return Interpolate1D(spectrum.XArray, spectrum.YArray, mzAxisNew); 
		}
		public double[] Interpolate1D(double[] x, double[] y, double[] xNew)
		{
			double[] yNew = new double[xNew.Length];
			double dx, dy, m, b;
			
			int xMaxIndex = x.Length - 1;
			int yMaxindex = y.Length - 1;
			int xNewSize = xNew.Length; 

			for(int i = 0; i < xNew.Length; i++)
			{
				int index = Array.BinarySearch(x, xNew[i]);
				if (index < 0) index = ~index;

				if (index > x.Length - 1) continue; 

				if(x[index] > xNew[i])
				{
					dx = index > 0 ? (x[index] - x[index - 1]) : (x[index + 1] - x[index]);
					dy = index > 0 ? (y[index] - y[index - 1]) : (y[index + 1] - y[index]);
				}else
				{
					dx = index < xMaxIndex ? (x[index + 1] - x[index]) : (x[index] - x[index - 1]);
					dy = index < xMaxIndex ? (y[index + 1] - y[index]) : (y[index] - y[index - 1]);
				}
				m = dy / dx;
				b = y[index] - x[index] * m;
				yNew[i] = xNew[i] * m + b; 
			}
			return yNew; 
		}
		public void NormalizeIntensity(ref double[] yarray)
		{
			double maxIntensity = yarray.Max();
			yarray = yarray.AsEnumerable().Select(i => i / maxIntensity).ToArray(); 
		}
		public double[,] PointSmoothing(double[,] matrix, int width)
		{
			for(int i = 0; i < matrix.GetLength(1); i++)
			{
				for(int j = 0; j < matrix.GetLength(0); j++)
				{
					int low = i - width;
					if (low < 0) low = 0;
					int high = i + width + 1;

					if (high > matrix.GetLength(1)) high = matrix.GetLength(1);

					double sum = 0;
					for (int k = low; k < high; k++)
					{
						sum += matrix[j, i] = sum / (double)(1 + 2 * width); 
					}
				}
			}
			return matrix; 
		}

	}

}
