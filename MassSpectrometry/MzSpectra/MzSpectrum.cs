// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MzSpectrum.cs) is part of MassSpectrometry.
//
// MassSpectrometry is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MzLibUtil;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MassSpectrometry
{
    public abstract class MzSpectrum<TPeak> : Spectrum<TPeak>, IMzSpectrum<TPeak>
        where TPeak : IMzPeak
    {

        #region Private Fields

        private static readonly List<List<double>> allMasses = new List<List<double>>();
        private static readonly List<List<double>> allIntensities = new List<List<double>>();

        private static readonly double[] mostIntenseMasses = new double[271];
        private static readonly int[] indicesOfMostIntense = new int[271];

        #endregion Private Fields

        #region Public Constructors

        static MzSpectrum()
        {
            // AVERAGINE
            const double averageC = 4.9384;
            const double averageH = 7.7583;
            const double averageO = 1.4773;
            const double averageN = 1.3577;
            const double averageS = 0.0417;

            int i = 0;
            do
            {
                i++;
                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * i));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * i));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * i));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * i));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * i));

                double fineRes = 0.5;
                var ye = IsotopicDistribution.GetDistribution(chemicalFormula, fineRes, 0);
                int firstCount = ye.Intensities.Count();
                do
                {
                    fineRes /= 2;
                    ye = IsotopicDistribution.GetDistribution(chemicalFormula, fineRes, 0);
                } while (ye.Intensities.Count() <= firstCount);
                fineRes *= 2;
                ye = IsotopicDistribution.GetDistribution(chemicalFormula, fineRes, 0);

                var masses = ye.Masses.ToList();
                var intensities = ye.Intensities.ToList();

                indicesOfMostIntense[i - 1] = intensities.IndexOf(intensities.Max());
                mostIntenseMasses[i - 1] = masses[indicesOfMostIntense[i - 1]];
                allMasses.Add(masses);
                allIntensities.Add(intensities);
            } while (mostIntenseMasses[i - 1] < 30000);
        }

        #endregion Public Constructors

        #region Protected Constructors

        protected MzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        protected MzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Protected Constructors

        #region Public Properties

        new public MzRange Range
        {
            get
            {
                return new MzRange(FirstX, LastX);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static byte[] Get64Bitarray(IEnumerable<double> yArray)
        {
            var mem = new MemoryStream();
            foreach (var okk in yArray)
            {
                byte[] ok = BitConverter.GetBytes(okk);
                mem.Write(ok, 0, ok.Length);
            }
            mem.Position = 0;
            return mem.ToArray();
        }

        public byte[] Get64BitYarray()
        {
            return Get64Bitarray(YArray);
        }

        public byte[] Get64BitXarray()
        {
            return Get64Bitarray(XArray);
        }

        public void ReplaceXbyApplyingFunction(Func<IMzPeak, double> convertor)
        {
            for (int i = 0; i < Size; i++)
                XArray[i] = convertor(this[i]);
            peakWithHighestY = default(TPeak);
            peakList = new TPeak[Size];
        }

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Size);
        }

        // Mass tolerance must account for different isotope spacing!
        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, Func<IMzPeak, bool> peakFilter)
        {
            var isolatedMassesAndCharges = new List<Tuple<List<IMzPeak>, int, int, double, List<double>>>();

            foreach (var candidateForMostIntensePeak in Extract(theRange).Where(b => peakFilter(b)))
            {
                List<IMzPeak> bestListOfPeaks = new List<IMzPeak>();
                int bestChargeState = 1;
                int bestIndexOfMostIntense = 0;
                double bestMonoisotopicMass = 0;
                var bestListOfRatios = new List<double>();

                for (int chargeState = 1; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    var testMostIntenseMass = candidateForMostIntensePeak.Mz.ToMass(chargeState);

                    var indexOfMassInStaticArrays = Array.BinarySearch(mostIntenseMasses, testMostIntenseMass);
                    if (indexOfMassInStaticArrays < 0)
                        indexOfMassInStaticArrays = ~indexOfMassInStaticArrays;
                    if (indexOfMassInStaticArrays == mostIntenseMasses.Length)
                        indexOfMassInStaticArrays--;

                    var listOfPeaksForThisChargeState = new List<IMzPeak> { candidateForMostIntensePeak };
                    var listOfRatios = new List<double> { allIntensities[indexOfMassInStaticArrays][indicesOfMostIntense[indexOfMassInStaticArrays]] / candidateForMostIntensePeak.Intensity };
                    // Assuming the test peak is most intense...
                    // Try to find the rest of the isotopes!

                    int indexOfLowestFoundMass = indicesOfMostIntense[indexOfMassInStaticArrays];
                    int indexDown = indicesOfMostIntense[indexOfMassInStaticArrays] - 1;
                    int indexUp = indicesOfMostIntense[indexOfMassInStaticArrays] + 1;
                    while (indexDown >= 0 || indexUp < allIntensities[indexOfMassInStaticArrays].Count)
                    {
                        int indexToLookAt;
                        if (indexDown < 0)
                        {
                            indexToLookAt = indexUp;
                            indexUp++;
                        }
                        else if (indexUp == allIntensities[indexOfMassInStaticArrays].Count)
                        {
                            indexToLookAt = indexDown;
                            indexDown--;
                        }
                        else if (allIntensities[indexOfMassInStaticArrays][indexDown] > allIntensities[indexOfMassInStaticArrays][indexUp])
                        {
                            indexToLookAt = indexDown;
                            indexDown--;
                        }
                        else
                        {
                            indexToLookAt = indexUp;
                            indexUp++;
                        }

                        double theorMassThatTryingToFind = testMostIntenseMass + allMasses[indexOfMassInStaticArrays][indexToLookAt] - allMasses[indexOfMassInStaticArrays][indicesOfMostIntense[indexOfMassInStaticArrays]];
                        var closestPeakToTheorMass = GetClosestPeak(theorMassThatTryingToFind.ToMz(chargeState));
                        if (Math.Abs(closestPeakToTheorMass.Mz.ToMass(chargeState) - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 <= deconvolutionTolerancePpm
                            && Peak2satisfiesRatio(allIntensities[indexOfMassInStaticArrays][indicesOfMostIntense[indexOfMassInStaticArrays]], allIntensities[indexOfMassInStaticArrays][indexToLookAt], candidateForMostIntensePeak, closestPeakToTheorMass, intensityRatioLimit))
                        {
                            // Found a match to an isotope peak for this charge state!
                            listOfPeaksForThisChargeState.Add(closestPeakToTheorMass);
                            listOfRatios.Add(allIntensities[indexOfMassInStaticArrays][indexToLookAt] / closestPeakToTheorMass.Intensity);
                            indexOfLowestFoundMass = Math.Min(indexToLookAt, indexOfLowestFoundMass);
                        }
                        else
                            break;
                    }
                    if (listOfPeaksForThisChargeState.Count > bestListOfPeaks.Count)
                    {
                        bestListOfPeaks = listOfPeaksForThisChargeState;
                        bestChargeState = chargeState;
                        bestMonoisotopicMass = listOfPeaksForThisChargeState.Min(b => b.Mz).ToMass(chargeState) - allMasses[indexOfMassInStaticArrays][indexOfLowestFoundMass] + allMasses[indexOfMassInStaticArrays][0];
                        bestIndexOfMostIntense = indicesOfMostIntense[indexOfMassInStaticArrays];
                        bestListOfRatios = listOfRatios;
                    }
                }

                if (bestListOfPeaks.Count >= 2)
                    isolatedMassesAndCharges.Add(new Tuple<List<IMzPeak>, int, int, double, List<double>>(bestListOfPeaks, bestChargeState, bestIndexOfMostIntense, bestMonoisotopicMass, bestListOfRatios));
            }

            HashSet<double> seen = new HashSet<double>();
            foreach (var ok in isolatedMassesAndCharges.OrderByDescending(b => b.Item1.Count - MathNet.Numerics.Statistics.Statistics.StandardDeviation(b.Item5)))
            {
                if (seen.Overlaps(ok.Item1.Select(b => b.Mz)))
                    continue;
                foreach (var ah in ok.Item1.Select(b => b.Mz))
                    seen.Add(ah);
                yield return new IsotopicEnvelope(ok);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private bool Peak2satisfiesRatio(double peak1theorIntensity, double peak2theorIntensity, IMzPeak peak1, IMzPeak peak2, double intensityRatio)
        {
            var comparedShouldBe = peak1.Intensity / peak1theorIntensity * peak2theorIntensity;

            if (peak2.Intensity < comparedShouldBe / intensityRatio || peak2.Intensity > comparedShouldBe * intensityRatio)
                return false;

            return true;
        }

        #endregion Private Methods

    }
}