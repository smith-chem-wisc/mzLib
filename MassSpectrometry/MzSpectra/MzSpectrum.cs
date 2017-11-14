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

        private const int numAveraginesToGenerate = 1500;
        private static readonly double[][] allMasses = new double[numAveraginesToGenerate][];
        private static readonly double[][] allIntensities = new double[numAveraginesToGenerate][];
        private static readonly double[] mostIntenseMasses = new double[numAveraginesToGenerate];
        private static readonly double[] diffToMonoisotopic = new double[numAveraginesToGenerate];

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

            const double fineRes = 0.125;
            const double minRes = 1e-8;

            for (int i = 0; i < numAveraginesToGenerate; i++)
            {
                double averagineMultiplier = (i + 1) / 2.0;
                //Console.Write("numAveragines = " + numAveragines);
                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * averagineMultiplier));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * averagineMultiplier));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * averagineMultiplier));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * averagineMultiplier));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * averagineMultiplier));

                {
                    var chemicalFormulaReg = chemicalFormula;
                    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                    var masses = ye.Masses.ToArray();
                    var intensities = ye.Intensities.ToArray();
                    Array.Sort(intensities, masses);
                    Array.Reverse(intensities);
                    Array.Reverse(masses);

                    mostIntenseMasses[i] = masses[0];
                    diffToMonoisotopic[i] = masses[0] - chemicalFormulaReg.MonoisotopicMass;
                    allMasses[i] = masses;
                    allIntensities[i] = intensities;
                }
            }
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

        public static byte[] Get64Bitarray(IEnumerable<double> array)
        {
            var mem = new MemoryStream();
            foreach (var okk in array)
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

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Size);
        }

        // Mass tolerance must account for different isotope spacing!
        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit)
        {
            var isolatedMassesAndCharges = new List<IsotopicEnvelope>();

            foreach (var candidateForMostIntensePeak in ExtractIndices(theRange.Minimum, theRange.Maximum))
            {
                IsotopicEnvelope bestIsotopeEnvelopeForThisPeak = null;

                var candidateForMostIntensePeakMz = XArray[candidateForMostIntensePeak];
                //Console.WriteLine("candidateForMostIntensePeakMz: " + candidateForMostIntensePeakMz);
                var candidateForMostIntensePeakIntensity = YArray[candidateForMostIntensePeak];

                for (int chargeState = minAssumedChargeState; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    //Console.WriteLine(" chargeState: " + chargeState);
                    var testMostIntenseMass = candidateForMostIntensePeakMz.ToMass(chargeState);

                    var massIndex = Array.BinarySearch(mostIntenseMasses, testMostIntenseMass);
                    if (massIndex < 0)
                        massIndex = ~massIndex;
                    if (massIndex == mostIntenseMasses.Length)
                    {
                        //Console.WriteLine("Breaking  because mass is too high: " + testMostIntenseMass);
                        break;
                    }
                    //Console.WriteLine("  massIndex: " + massIndex);

                    var listOfPeaks = new List<(double, double)> { (candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity) };
                    var listOfRatios = new List<double> { allIntensities[massIndex][0] / candidateForMostIntensePeakIntensity };
                    // Assuming the test peak is most intense...
                    // Try to find the rest of the isotopes!

                    double differenceBetweenTheorAndActual = testMostIntenseMass - mostIntenseMasses[massIndex];
                    double totalIntensity = candidateForMostIntensePeakIntensity;
                    for (int indexToLookAt = 1; indexToLookAt < allIntensities[massIndex].Length; indexToLookAt++)
                    {
                        //Console.WriteLine("   indexToLookAt: " + indexToLookAt);
                        double theorMassThatTryingToFind = allMasses[massIndex][indexToLookAt] + differenceBetweenTheorAndActual;
                        //Console.WriteLine("   theorMassThatTryingToFind: " + theorMassThatTryingToFind);
                        //Console.WriteLine("   theorMassThatTryingToFind.ToMz(chargeState): " + theorMassThatTryingToFind.ToMz(chargeState));
                        var closestPeakToTheorMass = GetClosestPeakIndex(theorMassThatTryingToFind.ToMz(chargeState));
                        var closestPeakmz = XArray[closestPeakToTheorMass];
                        //Console.WriteLine("   closestPeakmz: " + closestPeakmz);
                        var closestPeakIntensity = YArray[closestPeakToTheorMass];
                        if (Math.Abs(closestPeakmz.ToMass(chargeState) - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 <= deconvolutionTolerancePpm
                            && Peak2satisfiesRatio(allIntensities[massIndex][0], allIntensities[massIndex][indexToLookAt], candidateForMostIntensePeakIntensity, closestPeakIntensity, intensityRatioLimit)
                            && !listOfPeaks.Contains((closestPeakmz, closestPeakIntensity)))
                        {
                            //Found a match to an isotope peak for this charge state!
                            //Console.WriteLine(" *   Found a match to an isotope peak for this charge state!");
                            //Console.WriteLine(" *   chargeState: " + chargeState);
                            //Console.WriteLine(" *   closestPeakmz: " + closestPeakmz);
                            listOfPeaks.Add((closestPeakmz, closestPeakIntensity));
                            totalIntensity += closestPeakIntensity;
                            listOfRatios.Add(allIntensities[massIndex][indexToLookAt] / closestPeakIntensity);
                        }
                        else
                            break;
                    }

                    var extrapolatedMonoisotopicMass = testMostIntenseMass - diffToMonoisotopic[massIndex]; // Optimized for proteoforms!!
                    var lowestMass = listOfPeaks.Min(b => b.Item1).ToMass(chargeState); // But may actually observe this small peak
                    var monoisotopicMass = Math.Abs(extrapolatedMonoisotopicMass - lowestMass) < 0.5 ? lowestMass : extrapolatedMonoisotopicMass;

                    IsotopicEnvelope test = new IsotopicEnvelope(listOfPeaks, monoisotopicMass, chargeState, totalIntensity, MathNet.Numerics.Statistics.Statistics.StandardDeviation(listOfRatios), massIndex);

                    if (listOfPeaks.Count >= 2 && ScoreIsotopeEnvelope(test) > ScoreIsotopeEnvelope(bestIsotopeEnvelopeForThisPeak))
                    {
                        //Console.WriteLine("Better charge state is " + test.charge);
                        //Console.WriteLine("peaks: " + string.Join(",", listOfPeaks.Select(b => b.Item1)));
                        bestIsotopeEnvelopeForThisPeak = test;
                    }
                }

                if (bestIsotopeEnvelopeForThisPeak != null && bestIsotopeEnvelopeForThisPeak.peaks.Count >= 2)
                    isolatedMassesAndCharges.Add(bestIsotopeEnvelopeForThisPeak);
            }

            HashSet<double> seen = new HashSet<double>();
            foreach (var ok in isolatedMassesAndCharges.OrderByDescending(b => ScoreIsotopeEnvelope(b)))
            {
                if (seen.Overlaps(ok.peaks.Select(b => b.mz)))
                    continue;
                foreach (var ah in ok.peaks.Select(b => b.mz))
                    seen.Add(ah);
                yield return ok;
            }
        }

        #endregion Public Methods

        #region Private Methods

        private double ScoreIsotopeEnvelope(IsotopicEnvelope b)
        {
            if (b == null)
                return 0;
            return b.totalIntensity / Math.Pow(b.stDev, 0.13) * Math.Pow(b.peaks.Count, 0.4) / Math.Pow(b.charge, 0.06);
        }

        private bool Peak2satisfiesRatio(double peak1theorIntensity, double peak2theorIntensity, double peak1intensity, double peak2intensity, double intensityRatio)
        {
            var comparedShouldBe = peak1intensity / peak1theorIntensity * peak2theorIntensity;

            if (peak2intensity < comparedShouldBe / intensityRatio || peak2intensity > comparedShouldBe * intensityRatio)
                return false;

            return true;
        }

        #endregion Private Methods
    }
}