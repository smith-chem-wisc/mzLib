using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InSourceDecay
{
	public class MatrixTuple
	{
		/*
		// Item1 = m/z_i; Item2 = m/z_j; Item3 = value at i,j
		public List<(double, double, double)> TuplizedMatrix { get; set; }
		public int FirstScan { get; set; }
		public int LastScan { get; set; }

		public MatrixTuple(double[,] matrix)
		{
			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				for(int j = 0; j < matrix.GetLength(1); j++)
				{
					Tuple<double, double, double> resultTuple = new Tuple<double, double, double>((double)i, (double)j, matrix[i, j]); 
				}
			}
		}
		public MatrixTuple(double[,] matrix, double[] matrixValues, bool filterZero = false)
		{
			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				for(int j = 0; j < matrix.GetLength(1); j++)
				{
					if (filterZero)
					{
						if (matrix[i, j] != (double)0)
						{
							(double, double, double) returnTuple = (matrixValues[i], matrixValues[j], matrix[i, j]);
							TuplizedMatrix.Add(returnTuple);
						}
						else
						{
							continue;
						}
					}
					else
					{
						(double, double, double) returnTuple = (matrixValues[i], matrixValues[j], matrix[i, j]);
						TuplizedMatrix.Add(returnTuple);
					}
				}
			}
		}
		public MatrixTuple(PartialCovarianceMatrix pCovMat, bool filterZero = false, bool stdevMatrix = false)
		{
			double[,] matrix; 
			
			if (stdevMatrix)
			{
				matrix = pCovMat.SignalToNoiseMatrix;
			}
			else
			{
				matrix = pCovMat.PCovMatrix; 
			}
			double[] matrixValues = pCovMat.IndexMZ;

			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				for (int j = 0; j < matrix.GetLength(1); j++)
				{
					if (filterZero)
					{
						if (matrix[i, j] != (double)0)
						{
							(double, double, double) returnTuple = (matrixValues[i], matrixValues[j], matrix[i, j]);
							TuplizedMatrix.Add(returnTuple);
						}
						else
						{
							continue;
						}
					}
					else
					{
						(double, double, double) returnTuple = (matrixValues[i], matrixValues[j], matrix[i, j]);
						TuplizedMatrix.Add(returnTuple);
					}
				}
			}
		}

		public MatrixTuple(AveragedScan matrix)
		{
			// TO DO: Add conversion from AveragedScan to MatrixTuple

		} 
		*/
	} 
}
