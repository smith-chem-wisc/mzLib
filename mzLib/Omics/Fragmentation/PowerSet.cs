using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.Fragmentation
{
    /// <summary>
    /// Provides methods to compute the power set of a given set of items.
    /// </summary>
    public class PowerSet
    {
        /// <summary>
        /// Return a list of the unique sums of all subsets of the input list of numbers, where the size of the subsets is less than or equal to the specified SubsetSize.
        /// </summary>
        /// <param name="numberList"> The input list of numbers to compute subset sums from.</param>
        /// <param name="subsetSize"> The maximum size of subsets to consider when computing sums. Only subsets with size less than or equal to this value will be included in the results.</param>
        /// <returns> A list of unique sums of subsets of the input list, where each subset's size is less than or equal to the specified SubsetSize.</returns>
        public static List<double> UniqueSubsetSums(IReadOnlyList<double> numberList, int subsetSize) 
        {
            if (subsetSize < 0) 
            {
                throw new ArgumentException("Subset size must be non-negative.");
            }

            switch (numberList.Count)
                {
                case 0:
                    return new List<double> { 0 }; // The sum of the empty set is defined as 0
                case > 15:
                    throw new ArgumentException("Input list is too large. Maximum supported size is 15 to avoid excessive memory usage.");
                default:
                    break;
            }

            int numberCount = numberList.Count;
            Span<double> sums = stackalloc double[1 << numberCount]; // To store the sums of subsets, with a maximum of 2^n subsets for n numbers
            Span<byte> sizes = stackalloc byte[1 << numberCount]; // To track the size of each subset

            int count = 1; // The index of the sums array
            sums[0] = 0; // Start with the empty set sum
            sizes[0] = 0; // Size of the empty set is 0

            for (int i = 0; i < numberCount; i++) 
            {
                double v = numberList[i];
                int currentCount = count;
                for (int j = 0; j < currentCount; j++) 
                {
                    int newSize = sizes[j] + 1;
                    if (newSize <= subsetSize) 
                    {
                        sums[count] = sums[j] + v;
                        sizes[count] = (byte)newSize;
                        count++;
                    }
                }
            }

            // Sort the sums, which will bring duplicates together for easy removal
            ArraySort(sums, count);

            //Remove duplicates by comparing adjacent sorted sums
            const double tolerance = 1e-9;
            List<double> result = new List<double>();
            for (int i = 0; i < count; i++) 
            {
                if (i == 0 || Math.Abs(sums[i] - sums[i - 1]) > tolerance) 
                {
                    result.Add(sums[i]);
                }
            }

            return result;
        }

        private static void ArraySort(Span<double> a, int n)
        {
            // Implement a sorting algorithm to sort the array in-place
            // This is a placeholder for the actual sorting implementation
            // insertion sort (fast for <=32 items)
            for (int i = 1; i < n; i++) // start at 1 since the first item is zero
            {
                double x = a[i];
                int j = i - 1;
                while (j >= 0 && a[j] > x)
                {
                    a[j + 1] = a[j];
                    j--;
                }
                a[j + 1] = x;
            }
        }


    }
}
