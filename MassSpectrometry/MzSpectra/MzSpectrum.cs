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

        private static readonly double[] mms = new double[] { 1.0029, 2.0052, 3.0077, 4.01, 5.012, 6.0139, 7.0154, 8.0164 };

        private static readonly List<Tuple<double, List<double>>> intensityFractions = new List<Tuple<double, List<double>>>();

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

            intensityFractions.Add(new Tuple<double, List<double>>(155, new List<double> { 0.915094568, 0.07782302, 0.006528797, 0.000289506 }));
            intensityFractions.Add(new Tuple<double, List<double>>(226, new List<double> { 0.88015657, 0.107467263, 0.011417303, 0.000730494 }));
            intensityFractions.Add(new Tuple<double, List<double>>(310, new List<double> { 0.837398069, 0.142430845, 0.01821746, 0.001683771 }));
            intensityFractions.Add(new Tuple<double, List<double>>(437, new List<double> { 0.777595132, 0.186958768, 0.031114269, 0.003704342, 0.000220493 }));
            intensityFractions.Add(new Tuple<double, List<double>>(620, new List<double> { 0.701235526, 0.238542629, 0.050903269, 0.008082801, 0.000985192 }));
            intensityFractions.Add(new Tuple<double, List<double>>(888, new List<double> { 0.602453248, 0.291899044, 0.084076553, 0.01790019, 0.002916629, 0.000410371 }));
            intensityFractions.Add(new Tuple<double, List<double>>(1243, new List<double> { 0.492328432, 0.333344333, 0.128351944, 0.035959923, 0.008063481, 0.001433271, 0.000195251 }));
            intensityFractions.Add(new Tuple<double, List<double>>(1797, new List<double> { 0.348495022, 0.336686099, 0.193731423, 0.082270917, 0.028068866, 0.008052644, 0.001907311, 0.000372359, 4.52281E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(2515, new List<double> { 0.229964408, 0.313975523, 0.238643189, 0.130654102, 0.056881604, 0.020732138, 0.006490044, 0.001706308, 0.000373761, 4.55951E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(3532, new List<double> { 0.12863395, 0.247015676, 0.254100853, 0.184302695, 0.104989402, 0.049731171, 0.020279668, 0.007267861, 0.002300006, 0.000619357, 9.64322E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(5019, new List<double> { 0.053526677, 0.145402081, 0.208920636, 0.209809764, 0.164605485, 0.107024765, 0.059770563, 0.029447041, 0.012957473, 0.005127018, 0.001845335, 0.000572486, 0.000115904 }));
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
                if (seen.Overlaps(ok.peaks.Select(b => b.Item1)))
                    continue;
                foreach (var ah in ok.peaks.Select(b => b.Item1))
                    seen.Add(ah);
                yield return ok;
            }
        }

        public IEnumerable<Tuple<List<IMzPeak>, int>> DeconvoluteOld(MzRange theRange, int maxAssumedChargeState, Tolerance massTolerance, double intensityRatio)
        {
            var isolatedMassesAndCharges = new List<Tuple<List<IMzPeak>, int>>();

            foreach (var peak in Extract(theRange))
            {
                // Always assume the current peak is a monoisotopic peak!

                List<IMzPeak> bestListOfPeaks = new List<IMzPeak>();
                int bestChargeState = 1;
                for (int chargeState = 1; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    var listOfPeaksForThisChargeState = new List<IMzPeak> { peak };
                    var mMass = peak.Mz.ToMass(chargeState);
                    for (int mm = 1; mm <= mms.Length; mm++)
                    {
                        double diffToNextMmPeak = mms[mm - 1];
                        double theorMass = mMass + diffToNextMmPeak;
                        var closestpeak = GetPeak(GetClosestPeakIndex(theorMass.ToMz(chargeState)));
                        if (massTolerance.Within(closestpeak.Mz.ToMass(chargeState), theorMass) && SatisfiesRatios(mMass, mm, peak, closestpeak, intensityRatio))
                        {
                            // Found a match to an isotope peak for this charge state!
                            listOfPeaksForThisChargeState.Add(closestpeak);
                        }
                        else
                            break;
                    }
                    if (listOfPeaksForThisChargeState.Count >= bestListOfPeaks.Count)
                    {
                        bestListOfPeaks = listOfPeaksForThisChargeState;
                        bestChargeState = chargeState;
                    }
                }
                if (bestListOfPeaks.Count >= 2)
                    isolatedMassesAndCharges.Add(new Tuple<List<IMzPeak>, int>(bestListOfPeaks, bestChargeState));
            }

            List<double> seen = new List<double>();
            while (isolatedMassesAndCharges.Any())
            {
                // Pick longest
                var longest = isolatedMassesAndCharges.OrderByDescending(b => b.Item1.Count).First();
                yield return longest;
                isolatedMassesAndCharges.Remove(longest);
                isolatedMassesAndCharges.RemoveAll(b => b.Item1.Intersect(longest.Item1).Any());
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

        private bool SatisfiesRatios(double mMass, int mm, IMzPeak ye, IMzPeak closestpeak, double intensityRatio)
        {
            double bestDiff = double.MaxValue;
            List<double> bestFracList = null;
            for (int i = 0; i < intensityFractions.Count; i++)
            {
                var diff = Math.Abs(mMass - intensityFractions[i].Item1);
                if (diff < bestDiff)
                {
                    bestDiff = diff;
                    bestFracList = intensityFractions[i].Item2;
                }
            }
            if (bestFracList == null || bestFracList.Count <= mm)
                return false;

            var theMM = bestFracList[0];
            var theCompared = bestFracList[mm];

            var comparedShouldBe = ye.Intensity / theMM * theCompared;

            if (closestpeak.Intensity < comparedShouldBe / intensityRatio || closestpeak.Intensity > comparedShouldBe * intensityRatio)
                return false;

            return true;
        }

        #endregion Private Methods
    }
}