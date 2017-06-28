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

        private const int numJ = 3;
        private const int numAveragines = 550;
        private static readonly double[][][] allMasses = new double[numJ][][];
        private static readonly double[][][] allIntensities = new double[numJ][][];
        private static readonly double[][] mostIntenseMasses = new double[numJ][];
        private static readonly double[][] diffToMonoisotopic = new double[numJ][];

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

            for (int j = 0; j < numJ; j++)
            {
                allMasses[j] = new double[numAveragines][];
                allIntensities[j] = new double[numAveragines][];
                mostIntenseMasses[j] = new double[numAveragines];
                diffToMonoisotopic[j] = new double[numAveragines];
            }

            for (int i = 0; i < numAveragines; i++)
            {
                double numAveragines = (i + 1) / 2.0;
                //Console.Write("numAveragines = " + numAveragines);
                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * numAveragines));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * numAveragines));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * numAveragines));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * numAveragines));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * numAveragines));

                {
                    var chemicalFormulaReg = chemicalFormula;
                    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                    var masses = ye.Masses.ToArray();
                    var intensities = ye.Intensities.ToArray();
                    Array.Sort(intensities, masses);
                    Array.Reverse(intensities);
                    Array.Reverse(masses);

                    mostIntenseMasses[0][i] = masses[0];
                    diffToMonoisotopic[0][i] = masses[0] - chemicalFormulaReg.MonoisotopicMass;
                    allMasses[0][i] = masses;
                    allIntensities[0][i] = intensities;
                }

                // Light
                {
                    int numberOfLysines = (int)(0.0582 * numAveragines);
                    ChemicalFormula chemicalFormulaLight = new ChemicalFormula(chemicalFormula);
                    chemicalFormulaLight.Add(PeriodicTable.GetElement(6)[13], 6 * numberOfLysines);
                    chemicalFormulaLight.Add(PeriodicTable.GetElement(6), -6 * numberOfLysines);
                    chemicalFormulaLight.Add(PeriodicTable.GetElement(7)[15], 2 * numberOfLysines);
                    chemicalFormulaLight.Add(PeriodicTable.GetElement(7), -2 * numberOfLysines);
                    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaLight, fineRes, minRes);
                    var masses = ye.Masses.ToArray();
                    var intensities = ye.Intensities.ToArray();
                    Array.Sort(intensities, masses);
                    Array.Reverse(intensities);
                    Array.Reverse(masses);

                    mostIntenseMasses[1][i] = masses[0];
                    diffToMonoisotopic[1][i] = masses[0] - chemicalFormulaLight.MonoisotopicMass;
                    allMasses[1][i] = masses;
                    allIntensities[1][i] = intensities;
                }

                // Heavy
                {
                    int numberOfLysines = (int)(0.0582 * numAveragines);
                    ChemicalFormula chemicalFormulaHeavy = new ChemicalFormula(chemicalFormula);
                    chemicalFormulaHeavy.Add(PeriodicTable.GetElement(1)[2], 8 * numberOfLysines);
                    chemicalFormulaHeavy.Add(PeriodicTable.GetElement(1), -8 * numberOfLysines);
                    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaHeavy, fineRes, minRes);
                    var masses = ye.Masses.ToArray();
                    var intensities = ye.Intensities.ToArray();
                    Array.Sort(intensities, masses);
                    Array.Reverse(intensities);
                    Array.Reverse(masses);

                    mostIntenseMasses[2][i] = masses[0];
                    diffToMonoisotopic[2][i] = masses[0] - chemicalFormulaHeavy.MonoisotopicMass;
                    allMasses[2][i] = masses;
                    allIntensities[2][i] = intensities;
                }

                //{
                //    ChemicalFormula chemicalFormulaMg = new ChemicalFormula(chemicalFormula);
                //    chemicalFormulaMg.Add(PeriodicTable.GetElement(12), 1);
                //    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaMg, fineRes, minRes);
                //    var masses = ye.Masses.ToArray();
                //    var intensities = ye.Intensities.ToArray();
                //    Array.Sort(intensities, masses);
                //    Array.Reverse(intensities);
                //    Array.Reverse(masses);

                //    mostIntenseMasses[3][i] = masses[0];
                //    diffToMonoisotopic[3][i] = masses[0] - chemicalFormulaMg.MonoisotopicMass;
                //    allMasses[3][i] = masses;
                //    allIntensities[3][i] = intensities;
                //}

                //{
                //    ChemicalFormula chemicalFormulaS = new ChemicalFormula(chemicalFormula);
                //    chemicalFormulaS.Add(PeriodicTable.GetElement(16), 1);
                //    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaS, fineRes, minRes);
                //    var masses = ye.Masses.ToArray();
                //    var intensities = ye.Intensities.ToArray();
                //    Array.Sort(intensities, masses);
                //    Array.Reverse(intensities);
                //    Array.Reverse(masses);

                //    mostIntenseMasses[4][i] = masses[0];
                //    diffToMonoisotopic[4][i] = masses[0] - chemicalFormulaS.MonoisotopicMass;
                //    allMasses[4][i] = masses;
                //    allIntensities[4][i] = intensities;
                //}

                //{
                //    ChemicalFormula chemicalFormulaCa = new ChemicalFormula(chemicalFormula);
                //    chemicalFormulaCa.Add(PeriodicTable.GetElement(20), 1);
                //    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaCa, fineRes, minRes);
                //    var masses = ye.Masses.ToArray();
                //    var intensities = ye.Intensities.ToArray();
                //    Array.Sort(intensities, masses);
                //    Array.Reverse(intensities);
                //    Array.Reverse(masses);

                //    mostIntenseMasses[5][i] = masses[0];
                //    diffToMonoisotopic[5][i] = masses[0] - chemicalFormulaCa.MonoisotopicMass;
                //    allMasses[5][i] = masses;
                //    allIntensities[5][i] = intensities;
                //}

                //// Fe

                //// Console.WriteLine();
                ////  Console.WriteLine("Fe");
                //ChemicalFormula chemicalFormulaFe = new ChemicalFormula(chemicalFormula);
                //chemicalFormulaFe.Add(PeriodicTable.GetElement(26), 1);
                //ye = IsotopicDistribution.GetDistribution(chemicalFormulaFe, fineRes, 0);
                //masses = ye.Masses.ToList();
                //intensities = ye.Intensities.ToList();
                ////  Console.WriteLine("masses = " + string.Join(" ", masses.Select(b => b.ToString("G9")).Take(30)));
                ////  Console.WriteLine("intensities = " + string.Join(" ", intensities.Select(b => b.ToString("G9")).Take(30)));

                //// Zn

                //// Console.WriteLine();
                ////  Console.WriteLine("Zn");
                //ChemicalFormula chemicalFormulaZn = new ChemicalFormula(chemicalFormula);
                //chemicalFormulaZn.Add(PeriodicTable.GetElement(30), 1);
                //ye = IsotopicDistribution.GetDistribution(chemicalFormulaZn, fineRes, 0);
                //masses = ye.Masses.ToList();
                //intensities = ye.Intensities.ToList();
                ////Console.WriteLine("masses = " + string.Join(" ", masses.Select(b => b.ToString("G9")).Take(30)));
                ////Console.WriteLine("intensities = " + string.Join(" ", intensities.Select(b => b.ToString("G9")).Take(30)));

                //indicesOfMostIntense[i - 1] = intensities.IndexOf(intensities.Max());
                //mostIntenseMasses[i - 1] = masses[indicesOfMostIntense[i - 1]];
                //allMasses.Add(masses);
                //allIntensities.Add(intensities);
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
            var isolatedMassesAndCharges = new List<IsotopicEnvelope>();

            foreach (var candidateForMostIntensePeak in Extract(theRange).Where(b => peakFilter(b)))
            {
                List<IMzPeak> bestListOfPeaks = new List<IMzPeak>();
                int bestChargeState = 1;
                double bestMonoisotopicMass = 0;
                double bestTotalIntensity = 0;
                double bestStDev = 0;
                int bestMassIndex = -1;
                int bestJ = -1;

                for (int chargeState = 1; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    var testMostIntenseMass = candidateForMostIntensePeak.Mz.ToMass(chargeState);

                    for (int j = 1; j < numJ; j++)
                    {
                        var massIndex = Array.BinarySearch(mostIntenseMasses[j], testMostIntenseMass);
                        if (massIndex < 0)
                            massIndex = ~massIndex;
                        if (massIndex == mostIntenseMasses[j].Length)
                            massIndex--;

                        var listOfPeaks = new List<IMzPeak> { candidateForMostIntensePeak };
                        var listOfRatios = new List<double> { allIntensities[j][massIndex][0] / candidateForMostIntensePeak.Intensity };
                        // Assuming the test peak is most intense...
                        // Try to find the rest of the isotopes!

                        double differenceBetweenTheorAndActual = testMostIntenseMass - mostIntenseMasses[j][massIndex];
                        double totalIntensity = candidateForMostIntensePeak.Intensity;
                        for (int indexToLookAt = 1; indexToLookAt < allIntensities[j][massIndex].Length; indexToLookAt++)
                        {
                            double theorMassThatTryingToFind = allMasses[j][massIndex][indexToLookAt] + differenceBetweenTheorAndActual;
                            var closestPeakToTheorMass = GetClosestPeak(theorMassThatTryingToFind.ToMz(chargeState));
                            if (Math.Abs(closestPeakToTheorMass.Mz.ToMass(chargeState) - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 <= deconvolutionTolerancePpm
                                && Peak2satisfiesRatio(allIntensities[j][massIndex][0], allIntensities[j][massIndex][indexToLookAt], candidateForMostIntensePeak, closestPeakToTheorMass, intensityRatioLimit))
                            {
                                // Found a match to an isotope peak for this charge state!
                                listOfPeaks.Add(closestPeakToTheorMass);
                                totalIntensity += closestPeakToTheorMass.Intensity;
                                listOfRatios.Add(allIntensities[j][massIndex][indexToLookAt] / closestPeakToTheorMass.Intensity);
                            }
                            else
                                break;
                        }
                        if (((totalIntensity - bestTotalIntensity) / totalIntensity) > 1e-6
                            || ((totalIntensity - bestTotalIntensity) / totalIntensity) > -1e-6 && MathNet.Numerics.Statistics.Statistics.StandardDeviation(listOfRatios) < bestStDev)
                        {
                            bestListOfPeaks = listOfPeaks;
                            bestMonoisotopicMass = testMostIntenseMass - diffToMonoisotopic[j][massIndex];
                            bestChargeState = chargeState;
                            bestTotalIntensity = totalIntensity;
                            bestStDev = MathNet.Numerics.Statistics.Statistics.StandardDeviation(listOfRatios);
                            bestMassIndex = massIndex;
                            bestJ = j;
                        }
                    }
                }

                if (bestListOfPeaks.Count >= 2)
                    isolatedMassesAndCharges.Add(new IsotopicEnvelope(bestListOfPeaks, bestMonoisotopicMass, bestChargeState, bestTotalIntensity, bestStDev, bestMassIndex, bestJ));
            }

            HashSet<double> seen = new HashSet<double>();
            foreach (var ok in isolatedMassesAndCharges.OrderByDescending(b => b.totalIntensity - b.stDev))
            {
                if (seen.Overlaps(ok.peaks.Select(b => b.Mz)))
                    continue;
                foreach (var ah in ok.peaks.Select(b => b.Mz))
                    seen.Add(ah);
                yield return ok;
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