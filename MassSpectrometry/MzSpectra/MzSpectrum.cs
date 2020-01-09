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
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MassSpectrometry
{
    public class MzSpectrum
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

        /// <summary>
        /// Shortreed's idea for top-down ms1 deconvolution (Jan. 6, 2020)
        /// Useful for non-isotopically resolved envelopes
        /// deconvolute the whole spectrum by treating every peak as charge 1, 2, 3, etc and recording the masses
        /// do a histogram to find which masses have multiple charge states: those are the real species
        /// Could possibly avoid doing the whole spectrum by just targeting the area of interest
        /// </summary>
        /// <param name="theRange"></param>
        /// <param name="minAssumedChargeState"></param>
        /// <param name="maxAssumedChargeState"></param>
        /// <param name="deconvolutionTolerancePpm"></param>
        /// <param name="intensityRatioLimit"></param>
        /// <returns></returns>
        //public IEnumerable<IsotopicEnvelope> TopDownDeconvolution(int minAssumedChargeState, int maxAssumedChargeState)
        //{
        //    //cycle through all charge states

        //}

        // Mass tolerance must account for different isotope spacing!
        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit)
        {
            if (Size == 0)
            {
                yield break;
            }

            var isolatedMassesAndCharges = new List<IsotopicEnvelope>();

            (int start, int end) indexes = ExtractIndices(theRange.Minimum, theRange.Maximum);
            //find most intense peak in the range
            double maxIntensity = 0;
            for(int index = indexes.start; index<indexes.end; index++)
            {
                if(YArray[index]>maxIntensity)
                {
                    maxIntensity = YArray[index];
                }
            }
            //foreach peak
            for (int candidateForMostIntensePeakIndex = indexes.start; candidateForMostIntensePeakIndex < indexes.end; candidateForMostIntensePeakIndex++)
            {
                double candidateForMostIntensePeakIntensity = YArray[candidateForMostIntensePeakIndex];
                if (candidateForMostIntensePeakIntensity * 20 >= maxIntensity) //ignore peptides that are over 10 times less intense than the most intense peak in the range.
                {
                    IsotopicEnvelope bestIsotopeEnvelopeForThisPeak = null;

                    double candidateForMostIntensePeakMz = XArray[candidateForMostIntensePeakIndex];

                    for (int chargeState = minAssumedChargeState; chargeState <= maxAssumedChargeState; chargeState++)
                    {
                        double testMostIntenseMass = candidateForMostIntensePeakMz.ToMass(chargeState);

                        int massIndex = Array.BinarySearch(mostIntenseMasses, testMostIntenseMass);
                        if (massIndex < 0)
                        {
                            massIndex = ~massIndex;
                        }
                        if (massIndex == mostIntenseMasses.Length)
                        {
                            break;
                        }

                        IsotopicEnvelope test = FindIsotopicEnvelope(massIndex, candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity,
                            testMostIntenseMass, chargeState, deconvolutionTolerancePpm, intensityRatioLimit);
                        if (test.Peaks.Count >= 2 && (bestIsotopeEnvelopeForThisPeak==null || test.Score > bestIsotopeEnvelopeForThisPeak.Score))
                        {
                            if (test.Charge >= 5) //look for other charge states if we suspect there to be multiple charge states
                            {
                                ObserveAdjacentChargeStates(test, candidateForMostIntensePeakMz, massIndex, deconvolutionTolerancePpm, intensityRatioLimit, minAssumedChargeState, maxAssumedChargeState);
                            }

                            bestIsotopeEnvelopeForThisPeak = test;
                        }
                    }

                    if (bestIsotopeEnvelopeForThisPeak != null)
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

        public IsotopicEnvelope FindIsotopicEnvelope(int massIndex, double candidateForMostIntensePeakMz, double candidateForMostIntensePeakIntensity, double testMostIntenseMass, int chargeState, double deconvolutionTolerancePpm, double intensityRatioLimit)
        {
            double[] theoreticalMasses = allMasses[massIndex];
            double[] theoreticalIntensities = allIntensities[massIndex];
            //add "most intense peak"
            var listOfObservedPeaks = new List<(double, double)> { (candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity) };
            var listOfRatios = new List<double> { theoreticalIntensities[0] / candidateForMostIntensePeakIntensity }; // theoreticalIntensities and theoreticalMasses are sorted by intensity, so first is most intense
            // Assuming the test peak is most intense...
            // Try to find the rest of the isotopes!
            double differenceBetweenTheorAndActualMass = testMostIntenseMass - mostIntenseMasses[massIndex]; //mass difference actual-theoretical
            double totalIntensity = candidateForMostIntensePeakIntensity;

            for (int indexToLookAt = 1; indexToLookAt < theoreticalIntensities.Length; indexToLookAt++) //cycle through all theoretical peaks in this envelope from most intense to least intense
            {
                double theorMassThatTryingToFind = theoreticalMasses[indexToLookAt] + differenceBetweenTheorAndActualMass; //get the mass to look for
                var closestPeakToTheorMass = GetClosestPeakIndex(theorMassThatTryingToFind.ToMz(chargeState)); //find the experimental peak
                var closestPeakmz = XArray[closestPeakToTheorMass.Value];
                var closestPeakIntensity = YArray[closestPeakToTheorMass.Value];
                //if the peak is within the deconvolution tolerance, has the correct intensity, and hasn't already been observed
                if (Math.Abs(closestPeakmz.ToMass(chargeState) - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 <= deconvolutionTolerancePpm
                    && Peak2satisfiesRatio(theoreticalIntensities[0], theoreticalIntensities[indexToLookAt], candidateForMostIntensePeakIntensity, closestPeakIntensity, intensityRatioLimit)
                    && !listOfObservedPeaks.Contains((closestPeakmz, closestPeakIntensity)))
                {
                    //Found a match to an isotope peak for this charge state!
                    listOfObservedPeaks.Add((closestPeakmz, closestPeakIntensity)); //add to observed list
                    totalIntensity += closestPeakIntensity; //add intensity
                    listOfRatios.Add(theoreticalIntensities[indexToLookAt] / closestPeakIntensity); //add ratio
                }
                else
                {
                    break;
                }
            }

            var extrapolatedMonoisotopicMass = testMostIntenseMass - diffToMonoisotopic[massIndex]; // Optimized for proteoforms!!
            var lowestMass = listOfObservedPeaks.Min(b => b.Item1).ToMass(chargeState); // But may actually observe this small peak
            var monoisotopicMass = Math.Abs(extrapolatedMonoisotopicMass - lowestMass) < 0.5 ? lowestMass : extrapolatedMonoisotopicMass;
            //This is a horrible way to determine the monoisotopic mass that does not use other isotopes for mass accuracy

            return new IsotopicEnvelope(listOfObservedPeaks, monoisotopicMass, chargeState, totalIntensity, MathNet.Numerics.Statistics.Statistics.StandardDeviation(listOfRatios), massIndex);
        }

        public void ObserveAdjacentChargeStates(IsotopicEnvelope originalEnvelope, double mostIntensePeakMz, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit, double minChargeToLookFor, double maxChargeToLookFor)
        {
            //look for the higher and lower charge state to confirm this thing exists
            //look lower
            int originalZ = originalEnvelope.Charge;
            double mostAbundantNeutralIsotope = mostIntensePeakMz.ToMass(originalZ);

            for (int lowerZ = originalZ - 1; lowerZ >= minChargeToLookFor; lowerZ--)
            {
                if (!FindChargeStateOfMass(originalEnvelope, lowerZ, mostAbundantNeutralIsotope, massIndex, deconvolutionTolerancePpm, intensityRatioLimit))
                {
                    break;
                }
            }

            for (int higherZ = originalZ + 1; higherZ<=maxChargeToLookFor; higherZ++)
            {
                if (!FindChargeStateOfMass(originalEnvelope, higherZ, mostAbundantNeutralIsotope, massIndex, deconvolutionTolerancePpm, intensityRatioLimit))
                {
                    break;
                }
            }
        }

        private bool FindChargeStateOfMass(IsotopicEnvelope originalEnvelope, int zToInvestigate, double mostAbundantNeutralIsotopeToInvestigate, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit)
        {
            double mostAbundantIsotopeMzForThisZ = mostAbundantNeutralIsotopeToInvestigate.ToMz(zToInvestigate);
            //let's go find that peak!
            int observedPeakIndex = Array.BinarySearch(XArray, mostAbundantIsotopeMzForThisZ);
            if (observedPeakIndex < 0)
            {
                observedPeakIndex = ~observedPeakIndex;
            }
            if (observedPeakIndex == XArray.Length)
            {
                return false;
            }
            IsotopicEnvelope test = FindIsotopicEnvelope(massIndex, XArray[observedPeakIndex], YArray[observedPeakIndex], mostAbundantNeutralIsotopeToInvestigate,
                zToInvestigate, deconvolutionTolerancePpm, intensityRatioLimit);
            
            //Want add this isotope score to the previous one. We have more peaks that support this mass than just the original isotopic envelope
            if (test.Score != 0)
            {
                originalEnvelope.AddToIsotopeScore(test.Score);
                return true;
            }
            else
            {
                return false;
            }
        }

        public void XCorrPrePreprocessing(double scanRangeMinMz, double scanRangeMaxMz, double precursorMz, double precursorDiscardRange = 1.5, double discreteMassBin = 1.0005079, double minimumAllowedIntensityRatioToBasePeak = 0.05)
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
            FilteringParams secondFilter = new FilteringParams(null, minimumAllowedIntensityRatioToBasePeak, nominalWindowWidthDaltons, null, true, false, false);

            MsDataFile.WindowModeHelper(ref genericIntensityArray, ref genericMzArray, secondFilter, genericMzArray.Min(), genericMzArray.Max(), true);

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
                scaledIntensities[i] = Math.Max(0, (genericIntensityArray[i] - (scaleValue/ (double)Math.Max(1, denominator))));
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

            this.YArray = intensites.ToArray();
            this.XArray = masses.ToArray();
            Array.Sort(this.XArray, this.YArray);
            XcorrProcessed = true;
        }

        public (int start, int end) ExtractIndices(double minX, double maxX)
        {
            int ind = Array.BinarySearch(XArray, minX);
            if (ind < 0)
            {
                ind = ~ind;
            }
            if(ind<Size && XArray[ind]<=maxX)
            {
                int start = ind;
                while (ind < Size && XArray[ind] <= maxX)
                {
                    ind++;
                }
                ind--;
                return (start, ind);
            }
            else
            {
                return (1,0);
            }
        }

        public int? GetClosestPeakIndex(double x)
        {
            if (Size == 0)
            {
                return null;
            }
            int index = Array.BinarySearch(XArray, x);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index >= Size)
            {
                return index - 1;
            }
            if (index == 0)
            {
                return index;
            }

            if (x - XArray[index - 1] > XArray[index] - x)
            {
                return index;
            }
            return index - 1;
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
            return XArray[GetClosestPeakIndex(x).Value];
        }

        public int NumPeaksWithinRange(double minX, double maxX)
        {
            int startingIndex = Array.BinarySearch(XArray, minX);
            if (startingIndex < 0)
            {
                startingIndex = ~startingIndex;
            }
            if (startingIndex >= Size)
            {
                return 0;
            }
            int endIndex = Array.BinarySearch(XArray, maxX);
            if (endIndex < 0)
            {
                endIndex = ~endIndex;
            }
            if (endIndex == 0)
            {
                return 0;
            }

            return endIndex - startingIndex;
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
            int ind = Array.BinarySearch(XArray, minX);
            if (ind < 0)
            {
                ind = ~ind;
            }
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

        private bool Peak2satisfiesRatio(double peak1theorIntensity, double peak2theorIntensity, double peak1intensity, double peak2intensity, double intensityRatio)
        {
            var comparedShouldBe = peak1intensity / peak1theorIntensity * peak2theorIntensity;

            if (peak2intensity < comparedShouldBe / intensityRatio || peak2intensity > comparedShouldBe * intensityRatio)
            {
                return false;
            }
            return true;
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