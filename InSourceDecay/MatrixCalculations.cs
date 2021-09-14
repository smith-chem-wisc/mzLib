using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics; 

namespace InSourceDecay
{
	public static class MatrixCalculations
	{
		public static double[,] Add(double[,] matrix1, double[,] matrix2)
		{
			if (matrix1.GetLength(0) != matrix2.GetLength(0) || matrix1.GetLength(1) != matrix2.GetLength(1))
			{
				throw new InvalidOperationException("Matrices are of unequal dimensions");
			}
			double[,] result = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
			for (int i = 0; i < matrix1.GetLength(0); i++)
			{
				for (int j = 0; j < matrix2.GetLength(1); j++)
				{
					result[i, j] = matrix1[i, j] + matrix2[i, j];
				}
			}
			return result;
		}
		public static double[,] Subtract(double[,] matrix1, double[,] matrix2)
		{
			double[,] result = new double[matrix1.GetLength(0), matrix2.GetLength(1)];
			if (matrix1.GetLength(0) == matrix2.GetLength(0) && matrix1.GetLength(1) == matrix2.GetLength(1))
			{				
				for (int i = 0; i < matrix1.GetLength(0); i++)
				{
					for (int j = 0; j < matrix2.GetLength(1); j++)
					{
						result[i, j] = matrix1[i, j] - matrix2[i, j];
					}
				}
			}
			else
			{
				throw new InvalidOperationException("Matrices are of unequal dimensions");
			}
			
			return result;
		}
		public static double[,] Multiply(double[,] matrix1, double[,] matrix2)
		{
			int matrix2rows = matrix2.GetLength(0);
			int matrix2cols = matrix2.GetLength(1);
			int matrix1rows = matrix1.GetLength(0);
			int matrix1cols = matrix1.GetLength(1);

			double[,] result = new double[matrix1rows, matrix2cols];
			if (matrix1cols != matrix2rows)
			{
				throw new InvalidOperationException("Matrix 1 rows does not equal Matrix 2 columns");
			}

			for (int i = 0; i < matrix1rows; i++)
			{
				for (int j = 0; j < matrix2cols; j++)
				{
					for (int k = 0; k < matrix1cols; k++)
					{
						result[i, j] += matrix1[i, k] * matrix2[k, j];
					}
				}
			}
			return result;
		}
		public static double[,] Multiply(double[,] matrix1, double scalar)
		{
			double[,] result = new double[matrix1.GetLength(0), matrix1.GetLength(1)];
			for (int i = 0; i < matrix1.GetLength(0); i++)
			{
				for (int j = 0; j < matrix1.GetLength(1); j++)
				{
					result[i, j] = matrix1[i, j] * scalar;
				}
			}
			return result;
		}
		public static double[,] Multiply(double[] vector1, double[] vector2)
		{
			double[,] result = new double[vector1.Length, vector2.Length];
			for (int i = 0; i < vector1.Length; i++)
			{
				for (int j = 0; j < vector2.Length; j++)
				{
					result[i, j] += vector1[i] * vector2[j];
				}
			}
			return result;
		}
		public static double[,] Transpose(double[,] matrix1)
		{
			double[,] result = new double[matrix1.GetLength(1), matrix1.GetLength(0)];
			for (int i = 0; i < matrix1.GetLength(0); i++)
			{
				for (int j = 0; j < matrix1.GetLength(1); j++)
				{
					result[j, i] = matrix1[i, j];
				}
			}
			return result;
		}
		public static double[,] Covariance(double[,] matrix)
		{
			int rows = matrix.GetLength(0);
			double[,] covMatrix = new double[rows, rows]; 

			// calculate the deviation matrix: 
			// x = X - 1*1'*X/n
			// where x is the deviation matrix, X is the input matrix, 1 is a 1 x n matrix filled with ones, 
			// 1' is its transpose and n is the number of rows in X. 
			double[,] onesMatrix = FillMatrixOf1s(rows, rows);

			double[,] step1 = Multiply(onesMatrix, matrix);
			double scalar1 = (double)1 / (double)rows;
			double[,] step2 = Multiply(step1, scalar1);
			double[,] deviationMatrix = Subtract(matrix, step2);
			double[,] deviationMatrixTranspose = Transpose(deviationMatrix);

			// Calculate the covariance matrix, which is V = X * X' / (n-1)
			double scalar2 = (double)1 / ((double)rows - 1);
			covMatrix = Multiply(Multiply(deviationMatrix, deviationMatrixTranspose), scalar2);
			return covMatrix;
		}

		// Creates an identify matrix of size rows x cols. 
		// Should probably be moved to the matrix class. 
		public static double[,] FillMatrixOf1s(int rows, int cols)
		{
			double[,] result = new double[rows, cols];
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					result[i, j] = 1;
				}
			}
			return result;
		}
		public static double[] Covariance(double[,] matrix, double[] vector)
		{
			int matrixRows = matrix.GetLength(0);
			int matrixCols = matrix.GetLength(1);

			double[] covarianceResult = new double[matrixRows];
			for (int i = 0; i < matrixRows; i++)
			{
				double[] tempArray = new double[matrixCols];
				for (int j = 0; j < matrixCols; j++)
				{
					tempArray[j] = matrix[i, j];
				}
				covarianceResult[i] = Covariance(tempArray, vector);
			}
			return covarianceResult;

		}
		public static double Covariance(double[] array1, double[] array2)
		{
			double n = (double)1 / ((double)array1.Length - (double)1);

			double array1mean = array1.Average();
			double array2mean = array2.Average();

			double sum1 = array1.Zip(array2, (x1, y1) => (x1 - array1mean) * (y1 - array2mean)).Sum();
			double result = sum1 * n;
			if (result != double.PositiveInfinity && result != double.NegativeInfinity)
			{
				return result;
			}
			else
			{
				return 0;
			}
		}
		// Takes a List<double[,]> and calculates a covariance matrix for each 2D array in the List. \
		// May need to be moved to another class
		public static List<double[,]> CovarianceMatrixWithAveragedScans(List<double[,]> averagedScans)
		{
			List<double[,]> covarianceMatrices = new List<double[,]>();
			foreach (double[,] matrix in averagedScans)
			{
				int rows = matrix.GetLength(0);
				int cols = matrix.GetLength(1);
				// calculate the deviation matrix
				double[,] onesMatrix = FillMatrixOf1s(rows, rows);

				double[,] deviationMatrix = Subtract(matrix, Multiply(Multiply(onesMatrix, matrix), (double)1 / rows));
				double[,] covMatrix = Multiply(Multiply(deviationMatrix, Transpose(deviationMatrix)), (double)1 / (rows - 1));
				covarianceMatrices.Add(covMatrix);
			}
			return covarianceMatrices;
		}
		public static double[,] CalculatePartialCovarianceTerm(double[,] matrix, double[] randomVariable, double allObsRVVariance)
		{
			int matrixRows = matrix.GetLength(0);
			// Partial covariance term = COV(Feature X, TIC)COV(Feature Y, TIC)/COV(TIC, TIC)
			double[] covXTIC = new double[matrixRows];
			covXTIC = Covariance(matrix, randomVariable);

			double[,] covArrayWithTic = Multiply(covXTIC, covXTIC);
			double[,] partialCovTermMatrix = Multiply(covArrayWithTic, (double)1 / allObsRVVariance);
			return partialCovTermMatrix;
		}
		public static double[,] CalculatePartialCovarianceMatrix(double[,] matrix, double[] randomVariable, double allObsVVariance)
		{
			double[,] covarianceMatrix = Covariance(matrix);
			double[,] partialCovarianceTermMatrix = CalculatePartialCovarianceTerm(matrix, randomVariable, allObsVVariance);
			double[,] partialCovarianceMatrix = Subtract(covarianceMatrix, partialCovarianceTermMatrix);
			return partialCovarianceMatrix; 
		}
		public static double[,] CalculatePartialCovarianceMatrix(AveragedScan aveScan)
		{
			double[,] covarianceMatrix = Covariance(aveScan.DataMatrix);
			double[,] partialCovarianceTermMatrix = CalculatePartialCovarianceTerm(aveScan.DataMatrix, aveScan.TICValues, aveScan.AllScansTICVariance
				
				);
			double[,] result = Subtract(covarianceMatrix, partialCovarianceTermMatrix);
			return result; 
		}
		public static double[] RowMeans(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(0)]; 

			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				double[] rowVector = new double[matrix.GetLength(1)]; 
				for(int j = 0; j < matrix.GetLength(1); j++)
				{
					rowVector[j] = matrix[i, j]; 
				}
				result[i] = rowVector.Average(); 
			}
			return result; 
		}
		public static double[] ColMeans(double[,] matrix)
		{
			double[] result = new double[matrix.GetLength(1)]; 

			for(int i = 0; i < matrix.GetLength(1); i++)
			{
				double[] colVector = new double[matrix.GetLength(0)]; 
				for(int j = 0; j < matrix.GetLength(0); j++)
				{
					result[i] = colVector.Average(); 
				}
			}
			return result; 
		}
		public static double StandardDeviation(double[] vector)
		{
			double result = 0;
			double sum = 0;  ; 
			double mean = vector.Average();
			double N = vector.Length; 

			for(int i = 0; i < vector.Length; i++)
			{
				sum += Math.Pow((vector[i] - mean), (double)2); 
			}
			result = Math.Sqrt((sum / N));
			return result; 

		}
		public static double[] GetRow(double[,] matrix, int index)
		{
			double[] result = new double[matrix.GetLength(0)]; 
			for(int j = 0; j < matrix.GetLength(0); j++)
			{
				result[j] = matrix[index, j]; 
			}
			return result; 
		}
		public static double[] GetCol(double[,] matrix, int index)
		{
			double[] result = new double[matrix.GetLength(1)]; 
			for(int i = 0; i < matrix.GetLength(1); i++)
			{
				result[i] = matrix[i, index]; 
			}
			return result; 
		}
	}
}
