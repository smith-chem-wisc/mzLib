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
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace MassSpectrometry
{
    public class MzSpectrum : IEquatable<MzSpectrum>
    {
        private const int numAveraginesToGenerate = 1500;
        private static readonly double[][] allMasses = new double[numAveraginesToGenerate][];
        private static readonly double[][] allIntensities = new double[numAveraginesToGenerate][];
        private static readonly double[] mostIntenseMasses = new double[numAveraginesToGenerate];
        private static readonly double[] diffToMonoisotopic = new double[numAveraginesToGenerate];

        private MzPeak[] peakList;
        private int? indexOfpeakWithHighestY;
        private double? sumOfAllY;

        public double[] XArray { get; private set; }
        public double[] YArray { get; private set; }

        public bool XcorrProcessed { get; private set; }

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

        public MzSpectrum(double[,] mzintensities)
        {
            var count = mzintensities.GetLength(1);

            XArray = new double[count];
            YArray = new double[count];
            Buffer.BlockCopy(mzintensities, 0, XArray, 0, sizeof(double) * count);
            Buffer.BlockCopy(mzintensities, sizeof(double) * count, YArray, 0, sizeof(double) * count);
            peakList = new MzPeak[Size];
        }

        public MzSpectrum(double[] mz, double[] intensities, bool shouldCopy)
        {
            if (shouldCopy)
            {
                XArray = new double[mz.Length];
                YArray = new double[intensities.Length];
                Array.Copy(mz, XArray, mz.Length);
                Array.Copy(intensities, YArray, intensities.Length);
            }
            else
            {
                XArray = mz;
                YArray = intensities;
            }
            peakList = new MzPeak[Size];
        }

        public MzRange Range
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                return new MzRange(FirstX.Value, LastX.Value);
            }
        }

        public double? FirstX
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                return XArray[0];
            }
        }

        public double? LastX
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                return XArray[Size - 1];
            }
        }

        public int Size { get { return XArray.Length; } }

        public int? IndexOfPeakWithHighesetY
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                if (!indexOfpeakWithHighestY.HasValue)
                {
                    indexOfpeakWithHighestY = Array.IndexOf(YArray, YArray.Max());
                }
                return indexOfpeakWithHighestY.Value;
            }
        }

        public double? YofPeakWithHighestY
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                return YArray[IndexOfPeakWithHighesetY.Value];
            }
        }

        public double? XofPeakWithHighestY
        {
            get
            {
                if (Size == 0)
                {
                    return null;
                }
                return XArray[IndexOfPeakWithHighesetY.Value];
            }
        }

        public double SumOfAllY
        {
            get
            {
                if (!sumOfAllY.HasValue)
                {
                    sumOfAllY = YArray.Sum();
                }
                return sumOfAllY.Value;
            }
        }

        #region Obsolete Deconvolution Methods, now part of ClassicDeconv

        /// <summary>
        /// Shortreed's idea for top-down ms1 deconvolution (Jan. 6, 2020)
        /// Potential utility for non-isotopically resolved envelopes
        /// deconvolute the whole spectrum by treating every peak as charge 1, 2, 3, etc and recording the masses
        /// do a histogram to find which masses have multiple charge states: those with high overlap are the real species
        /// Could possibly avoid doing the whole spectrum by just targeting areas of interest
        /// </summary>
        //public IEnumerable<IsotopicEnvelope> TopDownDeconvolution(int minAssumedChargeState, int maxAssumedChargeState)
        //{
        //    //cycle through all charge states

        //}

        // Mass tolerance must account for different isotope spacing!
        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int minAssumedChargeState,
            int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit)
        {
            //if no peaks, stop
            if (Size == 0)
            {
                yield break;
            }

            var isolatedMassesAndCharges = new List<IsotopicEnvelope>();

            (int start, int end) indexes = ExtractIndices(theRange.Minimum, theRange.Maximum);

            //find the most intense peak in the range
            double maxIntensity = 0;
            for (int index = indexes.start; index < indexes.end; index++)
            {
                if (YArray[index] > maxIntensity)
                {
                    maxIntensity = YArray[index];
                }
            }

            //go through each peak in the selected range and assume it is the most intense peak of its isotopic envelope (if it's not, it will hopefully get a low score)
            //cycle through possible charge states and select the one that has the best score (fit) with the averagine model
            for (int candidateForMostIntensePeakIndex = indexes.start; candidateForMostIntensePeakIndex < indexes.end; candidateForMostIntensePeakIndex++)
            {
                double candidateForMostIntensePeakIntensity = YArray[candidateForMostIntensePeakIndex];
                if (candidateForMostIntensePeakIntensity * 100 >= maxIntensity) //ignore peptides that are over 100 times less intense than the most intense peak in the range (heuristic from Top-Down yeast)
                {
                    IsotopicEnvelope bestIsotopeEnvelopeForThisPeak = null;

                    double candidateForMostIntensePeakMz = XArray[candidateForMostIntensePeakIndex];

                    //Find what charge states this peak might be based on the spacing of nearby peaks (assumes isotopic resolution)
                    HashSet<int> allPossibleChargeStates = new HashSet<int>();
                    for (int i = candidateForMostIntensePeakIndex + 1; i < XArray.Length; i++) //look at peaks of higher m/z
                    {
                        double deltaMass = XArray[i] - candidateForMostIntensePeakMz;
                        if (deltaMass < 1.1) //if we're past a Th spacing, we're no longer looking at the closest isotope
                        {
                            //get the lower bound charge state
                            int charge = (int)Math.Floor(1 / deltaMass); //e.g. deltaMass = 0.4 Th, charge is now 2 (but might be 3)
                            if (charge >= minAssumedChargeState && charge <= maxAssumedChargeState)
                            {
                                allPossibleChargeStates.Add(charge);
                            }
                            //get the upper bound charge state
                            charge++;
                            if (charge >= minAssumedChargeState && charge <= maxAssumedChargeState)
                            {
                                allPossibleChargeStates.Add(charge);
                            }
                        }
                        else
                        {
                            break;
                        }
                    }

                    //investigate the putative charge states
                    foreach (int chargeState in allPossibleChargeStates)
                    {
                        //get the mass of this peak assuming it's the charge we're looking at
                        double testMostIntenseMass = candidateForMostIntensePeakMz.ToMass(chargeState);

                        //get the index of the theoretical isotopic envelope for an averagine model that's close in mass
                        int massIndex = mostIntenseMasses.GetClosestIndex(testMostIntenseMass);

                        //create a list for each isotopic peak from this envelope. This is used to fine tune the monoisotopic mass and is populated in "FindIsotopicEnvelope"
                        List<double> monoisotopicMassPredictions = new List<double>();

                        //Look for other isotopes using the assumed charge state
                        IsotopicEnvelope putativeIsotopicEnvelope = FindIsotopicEnvelope(massIndex, candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity,
                            testMostIntenseMass, chargeState, deconvolutionTolerancePpm, intensityRatioLimit, monoisotopicMassPredictions);

                        if (putativeIsotopicEnvelope.Peaks.Count >= 2) //if there are at least two isotopes
                        {
                            //look for other charge states, using them for scoring and monoisotopic mass estimates
                            //need to use this method before comparing scores because it changes the score of the test envelope
                            int numOtherChargeStatesObserved = ObserveAdjacentChargeStates(putativeIsotopicEnvelope, candidateForMostIntensePeakMz, massIndex, deconvolutionTolerancePpm, intensityRatioLimit, minAssumedChargeState, maxAssumedChargeState, monoisotopicMassPredictions);

                            //is this the best charge state for this peak?
                            if ((bestIsotopeEnvelopeForThisPeak == null || putativeIsotopicEnvelope.Score > bestIsotopeEnvelopeForThisPeak.Score) && //and the score is better for this charge state than others
                             (putativeIsotopicEnvelope.Charge / 5 <= numOtherChargeStatesObserved)) //and if we suspect there to be multiple charge states and there are (higher the charge, more states expected, z=5, need 2 charge states, z=10, need 3 charge states, etc
                            {
                                putativeIsotopicEnvelope.SetMedianMonoisotopicMass(monoisotopicMassPredictions); //take the median mass from all of the isotopes (this is fine tuning!)
                                bestIsotopeEnvelopeForThisPeak = putativeIsotopicEnvelope;
                            }
                        }
                    }
                    if (bestIsotopeEnvelopeForThisPeak != null) //add this envelope (it might be wrong, but hopefully it has a low score and gets outscored later by the right thing)
                    {
                        isolatedMassesAndCharges.Add(bestIsotopeEnvelopeForThisPeak);
                    }
                }
            }

            HashSet<double> seen = new HashSet<double>();
            foreach (var ok in isolatedMassesAndCharges.OrderByDescending(b => b.Score))
            {
                if (seen.Overlaps(ok.Peaks.Select(b => b.mz)))
                {
                    continue;
                }
                foreach (var ah in ok.Peaks.Select(b => b.mz))
                {
                    seen.Add(ah);
                }
                yield return ok;
            }
        }

        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        public IsotopicEnvelope FindIsotopicEnvelope(int massIndex, double candidateForMostIntensePeakMz, double candidateForMostIntensePeakIntensity, double testMostIntenseMass, int chargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, List<double> monoisotopicMassPredictions)
        {
            double[] theoreticalMasses = allMasses[massIndex];
            double[] theoreticalIntensities = allIntensities[massIndex];
            //add "most intense peak"
            var listOfObservedPeaks = new List<(double, double)> { (candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity) };
            var listOfRatios = new List<double> { theoreticalIntensities[0] / candidateForMostIntensePeakIntensity }; // theoreticalIntensities and theoreticalMasses are sorted by intensity, so first is most intense
            // Assuming the test peak is most intense...
            // Try to find the rest of the isotopes!
            double differenceBetweenTheorAndActualMass = testMostIntenseMass - theoreticalMasses[0]; //mass difference actual-theoretical for the tallest peak (not necessarily the monoisotopic)
            double totalIntensity = candidateForMostIntensePeakIntensity;
            double monoisotopicMass = testMostIntenseMass - diffToMonoisotopic[massIndex]; //get the  monoisotopic by taking the most intense mass minus the expected mass difference between most intense and monoisotopic
            monoisotopicMassPredictions.Add(monoisotopicMass);
            for (int indexToLookAt = 1; indexToLookAt < theoreticalIntensities.Length; indexToLookAt++) //cycle through all theoretical peaks in this envelope from most intense to least intense
            {
                double theorMassThatTryingToFind = theoreticalMasses[indexToLookAt] + differenceBetweenTheorAndActualMass; //get the expected mass of the next most intense peak
                int closestPeakToTheorMass = GetClosestPeakIndex(theorMassThatTryingToFind.ToMz(chargeState)); //find the experimental peak for that mass
                double closestPeakmz = XArray[closestPeakToTheorMass];
                double closestPeakIntensity = YArray[closestPeakToTheorMass];
                double closestPeakMass = closestPeakmz.ToMass(chargeState);
                //if the peak is within the deconvolution tolerance, has the correct intensity, and hasn't already been observed
                if (Math.Abs(closestPeakMass - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 <= deconvolutionTolerancePpm
                    && Peak2satisfiesRatio(theoreticalIntensities[0], theoreticalIntensities[indexToLookAt], candidateForMostIntensePeakIntensity, closestPeakIntensity, intensityRatioLimit)
                    && !listOfObservedPeaks.Contains((closestPeakmz, closestPeakIntensity)))
                {
                    //Found a match to an isotope peak for this charge state!
                    listOfObservedPeaks.Add((closestPeakmz, closestPeakIntensity)); //add to observed list
                    totalIntensity += closestPeakIntensity; //add intensity
                    listOfRatios.Add(theoreticalIntensities[indexToLookAt] / closestPeakIntensity); //add ratio
                    double monoisotopicMassFromThisPeak = monoisotopicMass + closestPeakMass - theorMassThatTryingToFind;
                    monoisotopicMassPredictions.Add(monoisotopicMassFromThisPeak);
                }
                else
                {
                    break;
                }
            }

            return new IsotopicEnvelope(listOfObservedPeaks, monoisotopicMass, chargeState, totalIntensity, Statistics.StandardDeviation(listOfRatios), massIndex);
        }

        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        public int ObserveAdjacentChargeStates(IsotopicEnvelope originalEnvelope, double mostIntensePeakMz, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit, double minChargeToLookFor, double maxChargeToLookFor, List<double> monoisotopicMassPredictions)
        {
            //look for the higher and lower charge states using the proposed mass
            int numAdjacentChargeStatesObserved = 0;
            int originalZ = originalEnvelope.Charge;
            double mostAbundantNeutralIsotope = mostIntensePeakMz.ToMass(originalZ);

            //look at lower charge states until we don't see one or we hit the minimum
            for (int lowerZ = originalZ - 1; lowerZ >= minChargeToLookFor; lowerZ--)
            {
                if (FindChargeStateOfMass(originalEnvelope, lowerZ, mostAbundantNeutralIsotope, massIndex, deconvolutionTolerancePpm, intensityRatioLimit, monoisotopicMassPredictions))
                {
                    numAdjacentChargeStatesObserved++;
                }
                else
                {
                    break;
                }
            }

            //look at higher charge states until we don't see one or we hit the maximum
            for (int higherZ = originalZ + 1; higherZ <= maxChargeToLookFor; higherZ++)
            {
                if (FindChargeStateOfMass(originalEnvelope, higherZ, mostAbundantNeutralIsotope, massIndex, deconvolutionTolerancePpm, intensityRatioLimit, monoisotopicMassPredictions))
                {
                    numAdjacentChargeStatesObserved++;
                }
                else
                {
                    break;
                }
            }
            return numAdjacentChargeStatesObserved;
        }

        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        private bool FindChargeStateOfMass(IsotopicEnvelope originalEnvelope, int zToInvestigate, double mostAbundantNeutralIsotopeToInvestigate, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit, List<double> monoisotopicMassPredictions)
        {
            //we know the mass and the charge that we're looking for, just see if the expected m/z and its isotopes are there or not
            double mostAbundantIsotopeMzForThisZTheoretical = mostAbundantNeutralIsotopeToInvestigate.ToMz(zToInvestigate);
            //let's go find that peak!
            int observedPeakIndex = GetClosestPeakIndex(mostAbundantIsotopeMzForThisZTheoretical);

            double mostAbundantIsotopeMzObserved = XArray[observedPeakIndex];
            double mostAbundantIsotopeMassObserved = mostAbundantIsotopeMzObserved.ToMass(zToInvestigate);
            //make sure the expected and observed peak are within the mass tolerance
            if (Math.Abs(mostAbundantIsotopeMassObserved - mostAbundantNeutralIsotopeToInvestigate) / mostAbundantNeutralIsotopeToInvestigate * 1e6 <= deconvolutionTolerancePpm)
            {
                //get the isotopic envelope for this peak and add the masses from all the peaks of the envelope to the monoisotopic mass predictions
                IsotopicEnvelope test = FindIsotopicEnvelope(massIndex, mostAbundantIsotopeMzObserved, YArray[observedPeakIndex], mostAbundantIsotopeMassObserved,
                    zToInvestigate, deconvolutionTolerancePpm, intensityRatioLimit, monoisotopicMassPredictions);

                //Add this isotope score to the original charge state score. We now have more peaks that support this mass than just the original isotopic envelope
                //This is currently additive, which probably could be improved upon
                if (test.Score != 0)
                {
                    originalEnvelope.AggregateChargeStateScore(test);
                    return true;
                }
                else
                {
                    //remove the test mass from this failed isotopic envelope
                    monoisotopicMassPredictions.RemoveAt(monoisotopicMassPredictions.Count - 1);
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        public (int start, int end) ExtractIndices(double minX, double maxX)
        {
            int start = XArray.GetClosestIndex(minX, ArraySearchOption.Next);
            if (XArray[start] <= maxX)
            {
                int end = XArray.GetClosestIndex(maxX, ArraySearchOption.Previous);
                return (start, end);
            }
            else
            {
                return (1, 0);
            }
        }

        [Obsolete("Deconvolution Has been moved to the Deconvoluter Object")]
        private bool Peak2satisfiesRatio(double peak1theorIntensity, double peak2theorIntensity, double peak1intensity, double peak2intensity, double intensityRatio)
        {
            var comparedShouldBe = peak1intensity / peak1theorIntensity * peak2theorIntensity;

            if (peak2intensity < comparedShouldBe / intensityRatio || peak2intensity > comparedShouldBe * intensityRatio)
            {
                return false;
            }
            return true;
        }

        #endregion

        public static byte[] Get64Bitarray(IEnumerable<double> array)
        {
            using (var mem = new MemoryStream())
            {
                foreach (var okk in array)
                {
                    byte[] ok = BitConverter.GetBytes(okk);
                    mem.Write(ok, 0, ok.Length);
                }
                mem.Position = 0;
                var memory = mem.ToArray();
                return memory;
            }
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

        public void XCorrPrePreprocessing(double scanRangeMinMz, double scanRangeMaxMz,
            double precursorMz, double precursorDiscardRange = 1.5,
            double discreteMassBin = 1.0005079, double minimumAllowedIntensityRatioToBasePeak = 0.05)
        {
            //The discrete bin value 1.0005079 was from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            Array.Sort(XArray, YArray);
            int numberOfNominalMasses = (int)Math.Round((scanRangeMaxMz - scanRangeMinMz) / discreteMassBin, 0) + 1;

            double[] genericIntensityArray = new double[numberOfNominalMasses];
            double[] genericMzArray = new double[numberOfNominalMasses];

            double lowestMz = Math.Round(scanRangeMinMz / discreteMassBin, 0) * discreteMassBin;

            for (int i = 0; i < numberOfNominalMasses; i++)
            {
                genericMzArray[i] = i * discreteMassBin + lowestMz;
                genericIntensityArray[i] = 0;
            }

            int lowTheortetical = (int)Math.Round((precursorMz - precursorDiscardRange) / discreteMassBin, 0);
            int hiTheortetical = (int)Math.Round((precursorMz + precursorDiscardRange) / discreteMassBin, 0);

            //this leaves you with one intensity per nominal mass, even if they come in as more than one. Intensity is Square-rooted
            for (int i = 0; i < XArray.Length; i++)
            {
                //the nominalMass is really just the integer index
                int nominalMass = (int)Math.Round((XArray[i] - scanRangeMinMz) / discreteMassBin, 0);

                //this might have be be exclusive (i.e. get rid of the = sign) should eliminate unfragmented precursors
                if (nominalMass < numberOfNominalMasses && (XArray[i] <= lowTheortetical || XArray[i] >= hiTheortetical))
                {
                    if (!Double.IsNaN(Math.Sqrt(YArray[i])) && !Double.IsInfinity(Math.Sqrt(YArray[i])))
                    {
                        genericIntensityArray[nominalMass] = Math.Max(genericIntensityArray[nominalMass], Math.Sqrt(YArray[i]));
                    }
                }
            }

            //we've already filtered for when multiple mzs appear in a single nominal mass bin
            int nominalWindowWidthDaltons = (int)(Math.Round((scanRangeMaxMz - scanRangeMinMz) / 10d, 0));

            FilteringParams secondFilter = new FilteringParams(null,
                minimumAllowedIntensityRatioToBasePeak, nominalWindowWidthDaltons, null,
                true, false, false);

            WindowModeHelper.Run(ref genericIntensityArray, ref genericMzArray, secondFilter,
                genericMzArray.Min(), genericMzArray.Max(), true);

            Array.Sort(genericMzArray, genericIntensityArray);

            //Scale the intensities

            const int rangeEnd = 75; //from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            double[] scaledIntensities = new double[genericIntensityArray.Length];
            for (int i = 0; i < genericIntensityArray.Length; i++)
            {
                double scaleValue = 0;

                int low = Math.Max(0, i - rangeEnd);
                int high = Math.Min(genericIntensityArray.Length - 1, i + rangeEnd);
                int denominator = high - low + 1;

                for (int j = low; j <= high; j++)
                {
                    if (!double.IsNaN(genericIntensityArray[j]))
                    {
                        scaleValue += genericIntensityArray[j];
                    }
                    else
                    {
                        denominator--;
                    }
                }
                scaledIntensities[i] = Math.Max(0, (genericIntensityArray[i] - (scaleValue / (double)Math.Max(1, denominator))));
            }

            genericIntensityArray = scaledIntensities;

            List<double> intensites = new List<double>();
            List<double> masses = new List<double>();

            for (int i = 0; i < genericIntensityArray.Count(); i++)
            {
                if (genericIntensityArray[i] > 0.0001)
                {
                    intensites.Add(genericIntensityArray[i]);
                    masses.Add(genericMzArray[i]);
                }
            }

            YArray = intensites.ToArray();
            XArray = masses.ToArray();
            Array.Sort(this.XArray, this.YArray);
            XcorrProcessed = true;
        }

        public int GetClosestPeakIndex(double x)
        {
            return XArray.GetClosestIndex(x);
        }

        public void ReplaceXbyApplyingFunction(Func<MzPeak, double> convertor)
        {
            for (int i = 0; i < Size; i++)
            {
                XArray[i] = convertor(GetPeak(i));
            }
            peakList = new MzPeak[Size];
        }

        public virtual double[,] CopyTo2DArray()
        {
            double[,] data = new double[2, Size];
            const int size = sizeof(double);
            Buffer.BlockCopy(XArray, 0, data, 0, size * Size);
            Buffer.BlockCopy(YArray, 0, data, size * Size, size * Size);
            return data;
        }

        public double? GetClosestPeakXvalue(double x)
        {
            if (Size == 0)
            {
                return null;
            }
            return XArray.GetClosestValue(x);
        }

        /// <summary>
        /// Reports the number of peaks between minX and maxX, inclusive
        /// </summary>
        public int NumPeaksWithinRange(double minX, double maxX)
        {
            if (XArray.Last() < minX || XArray.First() > maxX)
                return 0;

            int startingIndex = XArray.GetClosestIndex(minX, ArraySearchOption.Next);
            int endIndex = XArray.GetClosestIndex(maxX, ArraySearchOption.Previous);
            return endIndex - startingIndex + 1;
        }

        public IEnumerable<MzPeak> FilterByNumberOfMostIntense(int topNPeaks)
        {
            var quantile = 1.0 - (double)topNPeaks / Size;
            quantile = Math.Max(0, quantile);
            quantile = Math.Min(1, quantile);
            double cutoffYvalue = YArray.Quantile(quantile);

            for (int i = 0; i < Size; i++)
            {
                if (YArray[i] >= cutoffYvalue)
                {
                    yield return GetPeak(i);
                }
            }
        }

        public IEnumerable<MzPeak> Extract(DoubleRange xRange)
        {
            return Extract(xRange.Minimum, xRange.Maximum);
        }

        public IEnumerable<MzPeak> Extract(double minX, double maxX)
        {
            int ind = XArray.GetClosestIndex(minX, ArraySearchOption.Next);
            while (ind < Size && XArray[ind] <= maxX)
            {
                yield return GetPeak(ind);
                ind++;
            }
        }

        public IEnumerable<MzPeak> FilterByY(double minY, double maxY)
        {
            for (int i = 0; i < Size; i++)
            {
                if (YArray[i] >= minY && YArray[i] <= maxY)
                {
                    yield return GetPeak(i);
                }
            }
        }

        public IEnumerable<MzPeak> FilterByY(DoubleRange yRange)
        {
            return FilterByY(yRange.Minimum, yRange.Maximum);
        }

        public double CalculateDotProductSimilarity(MzSpectrum spectrumToCompare, Tolerance tolerance)
        {
            //get arrays of m/zs and intensities
            double[] mz1 = XArray;
            double[] intensity1 = YArray;

            double[] mz2 = spectrumToCompare.XArray;
            double[] intensity2 = spectrumToCompare.YArray;

            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();
            int i = 0; //iterate through mz1
            int j = 0; //iterate through mz2

            //find where peaks match
            while (i != mz1.Length && j != mz2.Length)
            {
                double one = mz1[i];
                double two = mz2[j];
                if (tolerance.Within(one, two)) //if match
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(intensity2[j]);
                    i++;
                    j++;
                }
                else if (one > two)
                {
                    vector1.Add(0);
                    vector2.Add(intensity2[j]);
                    j++;
                }
                else //two>one
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(0);
                    i++;
                }
            }
            //wrap up leftover peaks
            for (; i < mz1.Length; i++)
            {
                vector1.Add(intensity1[i]);
                vector2.Add(0);
            }
            for (; j < mz2.Length; j++)
            {
                vector1.Add(0);
                vector2.Add(intensity2[j]);
            }

            //numerator of dot product
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator of dot product
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            //return dot product
            return numerator / denominator;
        }

        public bool Equals(MzSpectrum? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return (XArray.SequenceEqual(other.XArray) && YArray.SequenceEqual(other.YArray));
        }

        public override bool Equals(object? spectrumToCompare)
        {
            if (ReferenceEquals(null, spectrumToCompare)) return false;
            if (ReferenceEquals(this, spectrumToCompare)) return true;
            if (spectrumToCompare.GetType() != this.GetType()) return false;
            return Equals((MzSpectrum)spectrumToCompare);
        }

        public override int GetHashCode()
        {
            return ((IStructuralEquatable)XArray).GetHashCode(EqualityComparer<double>.Default) +
                   ((IStructuralEquatable)YArray).GetHashCode(EqualityComparer<double>.Default) +
                   XArray.Sum().GetHashCode() + YArray.Sum().GetHashCode();
        }


        private MzPeak GetPeak(int index)
        {
            if (peakList[index] == null)
            {
                peakList[index] = GeneratePeak(index);
            }
            return peakList[index];
        }

        private MzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index], YArray[index]);
        }
    }
}