using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry; 

namespace SimulatedData
{
    public class SimulatedProtein : SimulatedData
    {
        private IsotopicDistribution _isotopicDistribution { get; set; }
        private readonly record struct BinValue(int Bin, double Mz, double Intensity);

        public SimulatedProtein(IsotopicDistribution isotopicDistribution, 
            double mzLow, double mzHigh, int length, double spacing) 
            : base(length, mzLow, spacing)
        {
            _isotopicDistribution = isotopicDistribution;
            Xarray[0] = mzLow;
            for (int i = 1; i < Xarray.Length; i++)
            {
                Xarray[i] = spacing * i + mzLow; 
            }

            var arrays = ConsumeDistribution(_isotopicDistribution, spacing);
            Yarray = arrays.yVals.ToArray();
        }
        private class BinValueComparer : IComparer<BinValue>
        {
            public int Compare(BinValue x, BinValue y)
            {
                return x.Bin.CompareTo(y.Bin);
            }
        }
        private static List<BinValue> CreateBinValueList(double[] xArray, double[] yArray,
            double min, double binSize)
        {
            var binIndices = xArray
                .Select((w, i) =>
                    new { Index = i, Bin = (int)Math.Round((w - min) / binSize, MidpointRounding.AwayFromZero) });
            List<BinValue> binValues = new List<BinValue>();
            foreach (var bin in binIndices)
            {
                binValues.Add(new BinValue(bin.Bin, xArray[bin.Index], yArray[bin.Index]));
            }

            return binValues;
        }
        private (List<double> xVals, List<double> yVals) ConsumeDistribution(IsotopicDistribution distribution, double binSize)
        {
            // TODO: Fix the generation of the x-axis and the envelope positioning. 
            //double min = distribution.Masses.Min();
            double min = 500; 

            int numberOfBins = Length;
            List<BinValue> listBinValues = CreateBinValueList(distribution.Masses.ToArray(),
                distribution.Intensities.ToArray(), min, binSize);
            listBinValues.Sort((m, n) => m.Bin.CompareTo(n.Bin));

            List<double> xVals = new();
            List<double> yVals = new();

            for (int i = 0; i < numberOfBins; i++)
            {
                List<BinValue> binValRecord = new();
                int index = listBinValues.BinarySearch(new BinValue() { Bin = i }, new BinValueComparer());
                if (index < 0)
                {
                    index = ~index;
                }

                int k = 0;
                while (k + index <= listBinValues.Count - 1)
                {
                    // binary search gets just the first index that it finds, so you 
                    // need to check on the left and right side for elements that 
                    // match the bin value. 

                    if (index + k < listBinValues.Count && listBinValues[index + k].Bin == i)
                    {
                        binValRecord.Add(listBinValues[index + k]);
                        if (k == 0)
                        {
                            k++;
                            continue;
                        }
                    }

                    if (index - k > 0 && listBinValues[index - k].Bin == i)
                    {
                        binValRecord.Add(listBinValues[index - k]);
                    }

                    k++;
                }

                if (binValRecord.Count() > 1)
                {
                    xVals.Add(binValRecord.Average(m => m.Mz));
                    yVals.Add(binValRecord.Average(m => m.Intensity));
                    continue;
                }

                if (binValRecord.Count == 0)
                {
                    // TODO: x array values aren't handled correctly. Should be the value of the m/z axis at the particular bin, not a zero. 
                    xVals.Add(0);
                    yVals.Add(0);
                    continue;
                }

                xVals.Add(binValRecord.First().Mz);
                yVals.Add(binValRecord.First().Intensity);
            }

            return (xVals, yVals);
        }
        
    }
}
