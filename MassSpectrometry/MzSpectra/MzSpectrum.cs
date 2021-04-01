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
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class MzSpectrum
    {
        private static readonly double[][] allMasses;
        private static readonly double[][] allIntensities;
        private static readonly double[] mostIntenseMasses;
        private static readonly double[] diffToMonoisotopic;

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

            double maxMass = 100000;
            double averagineDaltonResolution = 50.0;

            double averagineMass = PeriodicTable.GetElement("C").PrincipalIsotope.AtomicMass * averageC
                + PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * averageH
                + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass * averageO
                + PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass * averageN
                + PeriodicTable.GetElement("S").PrincipalIsotope.AtomicMass * averageS;

            int numAveragines = (int)Math.Ceiling(maxMass / averagineDaltonResolution) + 1;

            allMasses = new double[numAveragines][];
            allIntensities = new double[numAveragines][];
            mostIntenseMasses = new double[numAveragines];
            diffToMonoisotopic = new double[numAveragines];

            for (int i = 1; i < numAveragines; i++)
            {
                double mass = i * averagineDaltonResolution;

                double averagines = mass / averagineMass;

                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * averagines));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * averagines));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * averagines));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * averagines));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * averagines));

                var chemicalFormulaReg = chemicalFormula;
                IsotopicDistribution isotopicDistribution = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                var masses = isotopicDistribution.Masses.ToArray();
                var intensities = isotopicDistribution.Intensities.ToArray();
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

        public IEnumerable<IsotopicEnvelope> Deconvolute(MzRange theRange, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm,
            double intensityRatioLimit = 3, int minPeaks = 2, double minCorrelationToAveragine = 0.4, double minFractionIntensityObserved = 0.4, double signalToNoiseRequired = 2.0,
            int scanNum = 0, List<string> testStuff = null)
        {
            // if no peaks in the scan, stop
            if (Size == 0)
            {
                yield break;
            }

            double noiseRange = 10;

            var tolerance = new PpmTolerance(deconvolutionTolerancePpm);
            List<IsotopicEnvelope> envelopeCandidates = new List<IsotopicEnvelope>();
            List<DeconvolutedPeak> deconvolutedPeaks = new List<DeconvolutedPeak>();
            HashSet<double> mzsClaimed = new HashSet<double>();
            HashSet<int> potentialChargeStates = new HashSet<int>();
            Dictionary<double, List<IsotopicEnvelope>> potentialEnvsPerPeak = new Dictionary<double, List<IsotopicEnvelope>>();

            for (int p = 0; p < this.XArray.Length; p++)
            {
                double mz = XArray[p];
                double intensity = YArray[p];

                // check to see if this peak is in the m/z deconvolution range
                if (mz < theRange.Minimum - noiseRange)
                {
                    continue;
                }
                else if (mz > theRange.Maximum + noiseRange)
                {
                    break;
                }

                // get rough list of charge states to check for based on m/z peaks around this peak
                potentialChargeStates.Clear();
                for (int i = p + 1; i < this.XArray.Length; i++)
                {
                    double potentialIsotopeMz = this.XArray[i];

                    if (potentialIsotopeMz > mz + 1.5)
                    {
                        break;
                    }

                    int z = Convert.ToInt32(1.0 / (potentialIsotopeMz - mz));

                    if (z >= minAssumedChargeState && z <= maxAssumedChargeState)
                    {
                        potentialChargeStates.Add(z);
                    }
                }

                // examine different charge state possibilities
                foreach (int z in potentialChargeStates)
                {
                    int massIndex = GetMassIndex(mz.ToMass(z));

                    if (massIndex >= mostIntenseMasses.Length)
                    {
                        continue;
                    }

                    IsotopicEnvelope deconEnv = GetIsotopicEnvelope(mz, intensity, z, tolerance, intensityRatioLimit, deconvolutedPeaks,
                        mzsClaimed, minFractionIntensityObserved, minCorrelationToAveragine);

                    if (deconEnv != null)
                    {
                        envelopeCandidates.Add(deconEnv);

                        foreach (var peak in deconEnv.Peaks)
                        {
                            if (!potentialEnvsPerPeak.ContainsKey(peak.mz))
                            {
                                potentialEnvsPerPeak.Add(peak.mz, new List<IsotopicEnvelope>());
                            }

                            potentialEnvsPerPeak[peak.mz].Add(deconEnv);
                        }
                    }
                }
            }

            //TODO: scoring function
            var envelopesOrderedByScore = envelopeCandidates
                .OrderByDescending(p => p.Peaks.Count).ThenByDescending(p => p.PearsonCorrelation);
            List<IsotopicEnvelope> parsimoniousEnvelopes = new List<IsotopicEnvelope>();

            // greedy algorithm. pick the smallest set of envelopes that could explain the spectrum's peaks, given filtering criteria
            IsotopicEnvelope nextEnvelope = envelopesOrderedByScore.FirstOrDefault();
            while (nextEnvelope != null)
            {
                var items = potentialEnvsPerPeak[nextEnvelope.Peaks.First().mz];
                IsotopicEnvelope chosenEnvelope = nextEnvelope;

                // score/decide...
                if (items.Count > 1 && nextEnvelope.MonoisotopicMass > 3000)
                {
                    int maxPeaks = nextEnvelope.Peaks.Count;
                    double maxCorr = nextEnvelope.PearsonCorrelation;
                    items.RemoveAll(p => p.Peaks.Count <= maxPeaks - 3 || p.PearsonCorrelation < maxCorr - 0.2);

                    chosenEnvelope = items.OrderByDescending(p => p.MonoisotopicMass).First();
                }

                // add to envelope list
                var redecon = GetIsotopicEnvelope(chosenEnvelope.Peaks.First().mz, chosenEnvelope.Peaks.First().intensity, chosenEnvelope.Charge,
                    tolerance, intensityRatioLimit, deconvolutedPeaks, mzsClaimed, minFractionIntensityObserved, minCorrelationToAveragine, chosenEnvelope.MassIndex);

                parsimoniousEnvelopes.Add(redecon);

                foreach (var peak in redecon.Peaks)
                {
                    mzsClaimed.Add(peak.mz);
                }

                // re-evaluate all the other isotopic envelopes in the list that overlap with the just-added envelope.
                // newly invalidated envelopes will be set to null
                for (int i = 0; i < envelopeCandidates.Count; i++)
                {
                    var candidate = envelopeCandidates[i];

                    if (candidate.Peaks.Any(p => mzsClaimed.Contains(p.mz)))
                    {
                        IsotopicEnvelope redeconEnv = GetIsotopicEnvelope(candidate.Peaks.First().mz, candidate.Peaks.First().intensity, candidate.Charge,
                            tolerance, intensityRatioLimit, deconvolutedPeaks, mzsClaimed, minFractionIntensityObserved, minCorrelationToAveragine, candidate.MassIndex);

                        envelopeCandidates[i] = redeconEnv;
                    }
                }

                // remove all invalid isotopic envelopes
                envelopeCandidates.RemoveAll(p => p == null);


                foreach (var item in potentialEnvsPerPeak)
                {
                    item.Value.Clear();
                }
                foreach (var item in envelopeCandidates)
                {
                    foreach (var peak in item.Peaks)
                    {
                        potentialEnvsPerPeak[peak.mz].Add(item);
                    }
                }

                // get the next-best envelope
                nextEnvelope = envelopesOrderedByScore.FirstOrDefault();
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

            List<double> nearbyNoisePeaks = new List<double>();
            foreach (IsotopicEnvelope parsEnv in parsimoniousEnvelopes)
            {
                double mz = parsEnv.Peaks.First().mz;

                // calculate noise level
                nearbyNoisePeaks.Clear();

                foreach (var noisePeak in noisePeaks.Where(p => p.mz < mz + noiseRange && p.mz > mz - noiseRange))
                {
                    nearbyNoisePeaks.Add(noisePeak.intensity);
                }
                nearbyNoisePeaks.Add(parsEnv.Peaks.Min(p => p.intensity));

                // calculate signal to noise
                parsEnv.Noise = nearbyNoisePeaks.Median();
                parsEnv.SN = parsEnv.Peaks.Max(p => p.intensity) / parsEnv.Noise;

                // filter by S/N
                if (parsEnv.SN > signalToNoiseRequired)
                {
                    noiseFilteredEnvelopes.Add(parsEnv);
                }
            }

            // return deconvoluted envelopes
            foreach (IsotopicEnvelope envelope in noiseFilteredEnvelopes
                .Where(p => p.Peaks.First().mz > theRange.Minimum && p.Peaks.First().mz < theRange.Maximum))
            {
                yield return envelope;
            }
        }

        public IsotopicEnvelope GetIsotopicEnvelope(double mz, double intensity, int z, Tolerance tolerance, double intensityRatioLimit, List<DeconvolutedPeak> deconvolutedPeaks,
            HashSet<double> alreadyClaimedMzs, double fracIntensityNeeded = 0, double pearsonCorrNeeded = 0, int? massIndex = null)
        {
            if (alreadyClaimedMzs.Contains(mz))
            {
                return null;
            }

            deconvolutedPeaks.Clear();
            double mass = mz.ToMass(z);

            if (massIndex == null)
            {
                // get the index of the LOWEST MASS averagine envelope that's close in mass
                massIndex = GetMassIndex(mass);

                if (massIndex >= mostIntenseMasses.Length)
                {
                    return null;
                }
            }

            var averagineEnvelopeMasses = allMasses[massIndex.Value];
            var averagineEnvelopeIntensities = allIntensities[massIndex.Value];
            double monoMass = mass - diffToMonoisotopic[massIndex.Value];

            int indOfMostIntense = Array.IndexOf(averagineEnvelopeIntensities, 1);

            // 1 is to the right, -1 is to the left in the envelope
            int isotopeDirection = 1;
            for (int i = indOfMostIntense; i < averagineEnvelopeMasses.Length && i >= 0; i += isotopeDirection)
            {
                double isotopeMassShift = averagineEnvelopeMasses[i] - averagineEnvelopeMasses[indOfMostIntense];
                double isotopeMassToLookFor = mass + isotopeMassShift;
                double theoreticalIsotopeMz = isotopeMassToLookFor.ToMz(z);
                double theoreticalIsotopeIntensity = averagineEnvelopeIntensities[i] * intensity;

                var peakIndex = this.GetClosestPeakIndex(theoreticalIsotopeMz);

                //TODO: look for other peaks in the scan that could be this isotope that meet the m/z tolerance
                var isotopeExperMz = XArray[peakIndex];
                var isotopeExperIntensity = YArray[peakIndex];

                double intensityRatio = isotopeExperIntensity / theoreticalIsotopeIntensity;

                if (tolerance.Within(isotopeExperMz.ToMass(z), isotopeMassToLookFor) // check mass tolerance
                    && intensityRatio < intensityRatioLimit && intensityRatio > 1 / intensityRatioLimit // check intensity tolerance
                    && !alreadyClaimedMzs.Contains(isotopeExperMz) // check to see if this peak has already been claimed by another envelope
                    && !deconvolutedPeaks.Select(p => p.ExperimentalMz).Contains(isotopeExperMz)) // check to see if this peak has already been used in this envelope
                {
                    deconvolutedPeaks.Add(new DeconvolutedPeak(isotopeExperMz, theoreticalIsotopeMz, z, isotopeExperIntensity,
                        theoreticalIsotopeIntensity, i, averagineEnvelopeIntensities[i]));
                }
                else
                {
                    if (isotopeDirection == 1)
                    {
                        isotopeDirection = -1;
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

            double maxTheorIntensity = deconvolutedPeaks.First().TheoreticalIntensity;
            for (int i = 0; i < averagineEnvelopeIntensities.Length; i++)
            {
                double expectedIsotopeIntensity = averagineEnvelopeIntensities[i] * maxTheorIntensity;
                expectedTotalIntensity += expectedIsotopeIntensity;
            }

            double fracIntensityObserved = sumIntensity / expectedTotalIntensity;

            // calculate correlation to averagine
            double corr = Correlation.Pearson(deconvolutedPeaks.Select(p => p.ExperimentalIntensity), deconvolutedPeaks.Select(p => p.TheoreticalIntensity));

            // this is just to save memory, but quality filtering can happen outside of this method after the envelope has been returned, if desired
            if (corr < pearsonCorrNeeded || fracIntensityObserved < fracIntensityNeeded)
            {
                return null;
            }

            // create + return the isotopic envelope object
            return new IsotopicEnvelope(deconvolutedPeaks.Select(p => (p.ExperimentalMz, p.ExperimentalIntensity)).ToList(),
                monoMass, z, deconvolutedPeaks.Sum(p => p.ExperimentalIntensity), massIndex.Value, corr, fracIntensityObserved);
        }

        public int GetMassIndex(double mass)
        {
            int massIndex = Array.BinarySearch(mostIntenseMasses, mass);

            if (massIndex < 0)
            {
                massIndex = ~massIndex;
            }

            //// mass index is often off by 1, need to check to see if previous index is a better match
            //if (massIndex > 0)
            //{
            //    double prevMass = mostIntenseMasses[massIndex - 1];
            //    double thisIndexMass = mostIntenseMasses[massIndex];

            //    if (Math.Abs(mass - prevMass) < Math.Abs(mass - thisIndexMass))
            //    {
            //        massIndex -= 1;
            //    }
            //}
            //for (int m = massIndex; m >= 0; m--)
            //{
            //    if (Math.Abs(mostIntenseMasses[m] - mass) < 0.5)
            //    {
            //        massIndex = m;
            //    }
            //}

            return massIndex;
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