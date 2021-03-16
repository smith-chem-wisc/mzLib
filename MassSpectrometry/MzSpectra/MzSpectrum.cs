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
using System.Diagnostics;
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
                    //Array.Sort(intensities, masses);
                    //Array.Reverse(intensities);
                    //Array.Reverse(masses);
                    var mostIntense = intensities.Max();
                    int indOfMostIntensePeak = Array.IndexOf(intensities, mostIntense);

                    for (int j = 0; j < intensities.Length; j++)
                    {
                        intensities[j] /= mostIntense;
                    }

                    mostIntenseMasses[i] = masses[indOfMostIntensePeak];
                    diffToMonoisotopic[i] = masses[indOfMostIntensePeak] - chemicalFormulaReg.MonoisotopicMass;
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
        /// Potential utility for non-isotopically resolved envelopes
        /// deconvolute the whole spectrum by treating every peak as charge 1, 2, 3, etc and recording the masses
        /// do a histogram to find which masses have multiple charge states: those with high overlap are the real species
        /// Could possibly avoid doing the whole spectrum by just targeting areas of interest
        /// </summary>
        //public IEnumerable<IsotopicEnvelope> TopDownDeconvolution(int minAssumedChargeState, int maxAssumedChargeState)
        //{
        //    //cycle through all charge states
        //}

        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit,
            int minPeaks = 3, double minCorrelationToAveragine = 0.6, double minFractionIntensityObserved = 0.6, double signalToNoiseRequired = 3.0)
        {
            // if no peaks in the scan, stop
            if (Size == 0)
            {
                yield break;
            }

            var tolerance = new PpmTolerance(deconvolutionTolerancePpm);
            List<IsotopicEnvelope> envelopeCandidates = new List<IsotopicEnvelope>();
            List<DeconvolutedPeak> deconvolutedPeaks = new List<DeconvolutedPeak>();
            HashSet<double> mzsClaimed = new HashSet<double>();

            for (int p = 0; p < this.XArray.Length; p++)
            {
                double mz = XArray[p];
                double intensity = YArray[p];

                // check to see if this peak is in the m/z deconvolution range
                if (mz < theRange.Minimum - 5)
                {
                    continue;
                }
                else if (mz > theRange.Maximum + 5)
                {
                    break;
                }

                // examine different charge state possibilities
                for (int z = minAssumedChargeState; z <= maxAssumedChargeState; z++)
                {
                    var deconEnv = GetIsotopicEnvelope(mz, intensity, z, tolerance, intensityRatioLimit, deconvolutedPeaks, mzsClaimed);

                    if (deconEnv != null
                        && deconEnv.Peaks.Count >= minPeaks
                        && deconEnv.PearsonCorrelation > minCorrelationToAveragine
                        && deconEnv.FracIntensityObserved > minFractionIntensityObserved)
                    {
                        envelopeCandidates.Add(deconEnv);
                    }
                }
            }

            // pick the smallest list of envelopes that could explain the spectrum's peaks, given filtering criteria
            List<IsotopicEnvelope> parsimoniousEnvelopes = new List<IsotopicEnvelope>();

            IEnumerable<IsotopicEnvelope> envelopesOrderedByScore = envelopeCandidates.OrderByDescending(p => p.Peaks.Count).ThenByDescending(p => p.PearsonCorrelation);
            IsotopicEnvelope env = envelopesOrderedByScore.FirstOrDefault();

            while (env != null)
            {
                if (env.Peaks.Any(p => mzsClaimed.Contains(p.mz)))
                {
                    // this envelope has some peaks that overlap with the peaks that have been claimed.
                    // rescore the envelope and re-add the new envelope to the list of envelope candidates
                    envelopeCandidates.Remove(env);
                    IsotopicEnvelope redeconEnv = GetIsotopicEnvelope(env.Peaks.First().mz, env.Peaks.First().intensity, env.Charge,
                        tolerance, intensityRatioLimit, deconvolutedPeaks, mzsClaimed);

                    if (redeconEnv != null
                        && redeconEnv.Peaks.Count >= minPeaks
                        && redeconEnv.PearsonCorrelation > minCorrelationToAveragine
                        && redeconEnv.FracIntensityObserved > minFractionIntensityObserved)
                    {
                        envelopeCandidates.Add(redeconEnv);
                    }
                }
                else
                {
                    // the envelope has no peaks that have been claimed - add it to the list of parsimonious envelopes
                    parsimoniousEnvelopes.Add(env);
                    env.Score = env.PearsonCorrelation * env.Peaks.Count;

                    foreach (var peak in env.Peaks)
                    {
                        mzsClaimed.Add(peak.mz);
                    }
                }

                env = envelopesOrderedByScore.FirstOrDefault();
            }

            // noise filtering
            List<IsotopicEnvelope> noiseFilteredEnvelopes = new List<IsotopicEnvelope>();
            List<(double mz, double intensity)> noisePeaks = new List<(double mz, double intensity)>();
            for (int i = 0; i < XArray.Length; i++)
            {
                if (!mzsClaimed.Contains(XArray[i]))
                {
                    noisePeaks.Add((XArray[i], YArray[i]));
                }
            }

            foreach (IsotopicEnvelope parsEnv in parsimoniousEnvelopes)
            {
                double mz = parsEnv.Peaks.First().mz;
                var noise = noisePeaks.Where(p => p.mz < mz + 50 && p.mz > mz - 50).Select(p => p.intensity).Median(); //TODO: fix, this can crash (calling .Median() on empty list)

                double sn = parsEnv.Peaks.Sum(p => p.intensity) / noise;

                parsEnv.SN = sn;
                parsEnv.Noise = noise;

                if (parsEnv.Peaks.First().intensity / noise > signalToNoiseRequired)
                {
                    noiseFilteredEnvelopes.Add(parsEnv);
                }
            }

            // return deconvoluted envelopes
            foreach (var item in noiseFilteredEnvelopes)
            {
                yield return item;
            }
        }

        public IsotopicEnvelope GetIsotopicEnvelope(double mz, double intensity, int z, Tolerance tolerance, double intensityRatioLimit, List<DeconvolutedPeak> deconvolutedPeaks,
            HashSet<double> alreadyClaimedMzs)
        {
            deconvolutedPeaks.Clear();
            double mass = mz.ToMass(z);

            // get the index of the averagine envelope that's close in mass
            int massIndex = Array.BinarySearch(mostIntenseMasses, mass);
            if (massIndex < 0)
            {
                massIndex = ~massIndex;
            }
            if (massIndex == mostIntenseMasses.Length)
            {
                return null;
            }

            // mass index is often off by 1, need to check to see if previous index is a better match
            if (massIndex > 0)
            {
                double prevMass = mostIntenseMasses[massIndex - 1];
                double thisIndexMass = mostIntenseMasses[massIndex];

                if (Math.Abs(mass - prevMass) < Math.Abs(mass - thisIndexMass))
                {
                    massIndex -= 1;
                }
            }

            var averagineEnvelopeMasses = allMasses[massIndex];
            var averagineEnvelopeIntensities = allIntensities[massIndex];
            double monoMass = mass - diffToMonoisotopic[massIndex];

            var mostIntense = averagineEnvelopeIntensities.Max();
            int indOfMostIntense = Array.IndexOf(averagineEnvelopeIntensities, mostIntense);

            int direction = 1;
            for (int i = indOfMostIntense; i < averagineEnvelopeMasses.Length && i >= 0; i += direction)
            {
                double isotopeMassShift = averagineEnvelopeMasses[i] - averagineEnvelopeMasses[indOfMostIntense];
                double isotopeMassToLookFor = mass + isotopeMassShift;
                double theoreticalIsotopeMz = isotopeMassToLookFor.ToMz(z);
                double theoreticalIsotopeIntensity = averagineEnvelopeIntensities[i] * intensity;

                var peakIndex = this.GetClosestPeakIndex(theoreticalIsotopeMz);

                //TODO: look for other peaks in the scan that could be this isotope that meet the m/z tolerance
                var isotopeExpMz = XArray[peakIndex];
                var isotopeExpIntensity = YArray[peakIndex];

                double expectedIntensity = averagineEnvelopeIntensities[i] * intensity;
                double intensityRatio = isotopeExpIntensity / expectedIntensity;

                if (tolerance.Within(isotopeExpMz.ToMass(z), isotopeMassToLookFor) // check mass tolerance
                    && !deconvolutedPeaks.Select(p => p.ExperimentalMz).Contains(isotopeExpMz) // check to see if this peak has already been counted
                    && intensityRatio < intensityRatioLimit && intensityRatio > 1 / intensityRatioLimit // check intensity tolerance
                    && !alreadyClaimedMzs.Contains(isotopeExpMz)) // check to see if this peak has already been claimed by another envelope
                {
                    deconvolutedPeaks.Add(new DeconvolutedPeak(isotopeExpMz, theoreticalIsotopeMz, z, isotopeExpIntensity, theoreticalIsotopeIntensity, i, averagineEnvelopeIntensities[i]));
                }
                else
                {
                    if (direction == 1)
                    {
                        direction = -1;
                        i = indOfMostIntense;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            if (deconvolutedPeaks.Count < 2)
            {
                return null;
            }

            // calculate % intensity missing
            double sumIntensity = deconvolutedPeaks.Sum(p => Math.Min(p.ExperimentalIntensity, p.TheoreticalIntensity));
            double expectedTotalIntensity = 0;

            double maxIntensity = deconvolutedPeaks.Max(p => p.TheoreticalIntensity);
            for (int i = 0; i < averagineEnvelopeIntensities.Length; i++)
            {
                double expectedIsotopeIntensity = averagineEnvelopeIntensities[i] * maxIntensity;
                expectedTotalIntensity += expectedIsotopeIntensity;
            }

            double fracIntensityObserved = sumIntensity / expectedTotalIntensity;

            // check correlation to averagine
            double corr = Correlation.Pearson(deconvolutedPeaks.Select(p => p.ExperimentalIntensity), deconvolutedPeaks.Select(p => p.TheoreticalIntensity));

            double envScore = corr * deconvolutedPeaks.Count;// * fracIntensityObserved;

            IsotopicEnvelope env = new IsotopicEnvelope(deconvolutedPeaks.Select(p => (p.ExperimentalMz, p.ExperimentalIntensity)).ToList(),
                monoMass, z, deconvolutedPeaks.Sum(p => p.ExperimentalIntensity), massIndex, envScore);

            env.PearsonCorrelation = corr;
            env.FracIntensityObserved = fracIntensityObserved;

            return env;
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

        public (int start, int end) ExtractIndices(double minX, double maxX)
        {
            int ind = Array.BinarySearch(XArray, minX);
            if (ind < 0)
            {
                ind = ~ind;
            }
            if (ind < Size && XArray[ind] <= maxX)
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
                return (1, 0);
            }
        }

        public int GetClosestPeakIndex(double x)
        {
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
            return XArray[GetClosestPeakIndex(x)];
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