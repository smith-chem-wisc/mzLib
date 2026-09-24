using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace MassSpectrometry
{
    internal class ClassicDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        private MzSpectrum spectrum;

        internal ClassicDeconvolutionAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        /// <summary>
        /// Override to deconvolute the spectra using the Classic Deconvolution algorithm
        /// </summary>
        /// <param name="spectrumToDeconvolute">spectrum to deconvolute</param>
        /// <param name="range">Range of peaks to deconvolute</param>
        /// <returns></returns>
        protected internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrumToDeconvolute, MzRange range)
        {
            var deconParams = DeconvolutionParameters as ClassicDeconvolutionParameters ?? throw new MzLibException("Deconvolution params and algorithm do not match");
            spectrum = spectrumToDeconvolute;
            //if no peaks, stop
            if (spectrum.Size == 0)
            {
                yield break;
            }

            var isolatedMassesAndCharges = new List<IsotopicEnvelope>();

            (int start, int end) indexes = ExtractIndices(range.Minimum, range.Maximum);

            //find the most intense peak in the range
            double maxIntensity = 0;
            for (int index = indexes.start; index < indexes.end; index++)
            {
                if (spectrum.YArray[index] > maxIntensity)
                {
                    maxIntensity = spectrum.YArray[index];
                }
            }

            //go through each peak in the selected range and assume it is the most intense peak of its isotopic envelope (if it's not, it will hopefully get a low score)
            //cycle through possible charge states and select the one that has the best score (fit) with the averagine model
            for (int candidateForMostIntensePeakIndex = indexes.start;
                 candidateForMostIntensePeakIndex < indexes.end;
                 candidateForMostIntensePeakIndex++)
            {
                double candidateForMostIntensePeakIntensity = spectrum.YArray[candidateForMostIntensePeakIndex];
                if (candidateForMostIntensePeakIntensity * 100 >=
                    maxIntensity) //ignore peptides that are over 100 times less intense than the most intense peak in the range (heuristic from Top-Down yeast)
                {
                    IsotopicEnvelope bestIsotopeEnvelopeForThisPeak = null;

                    double candidateForMostIntensePeakMz = spectrum.XArray[candidateForMostIntensePeakIndex];

                    //Find what charge states this peak might be based on the spacing of nearby peaks (assumes isotopic resolution)
                    HashSet<int> allPossibleChargeStates = new HashSet<int>();
                    for (int i = candidateForMostIntensePeakIndex + 1;
                         i < spectrum.XArray.Length;
                         i++) //look at peaks of higher m/z
                    {
                        double deltaMass = spectrum.XArray[i] - candidateForMostIntensePeakMz;
                        if (deltaMass <
                            1.1) //if we're past a Th spacing, we're no longer looking at the closest isotope
                        {
                            //get the lower bound charge state
                            int charge = 0;
                            // A duplicated m/z makes deltaMass exactly zero, so 1/deltaMass is infinite.
                            // The int cast of that changed with .NET 9 - it used to wrap to int.MinValue
                            // and now saturates to int.MaxValue - so guard rather than rely on either.
                            if (deltaMass == 0)
                            {
                                continue;
                            }

                            if (deconParams.Polarity == Polarity.Negative)
                            {
                                charge = (int)Math.Floor(-1 / deltaMass); //e.g. deltaMass = 0.4 Th, charge is now 2 (but might be 3)
                            }
                            else
                            {
                                charge = (int)Math.Floor(1 / deltaMass); //e.g. deltaMass = 0.4 Th, charge is now 2 (but might be 3)
                            }

                            if (charge >= deconParams.MinAssumedChargeState && charge <= deconParams.MaxAssumedChargeState)
                            {
                                allPossibleChargeStates.Add(charge);
                            }

                            //get the upper bound charge state
                            charge++;
                            if (charge >= deconParams.MinAssumedChargeState && charge <= deconParams.MaxAssumedChargeState)
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
                        int massIndex = AverageResidueModel.GetMostIntenseMassIndex(testMostIntenseMass);

                        //create a list for each isotopic peak from this envelope. This is used to fine tune the monoisotopic mass and is populated in "FindIsotopicEnvelope"
                        List<double> monoisotopicMassPredictions = new List<double>();

                        //Look for other isotopes using the assumed charge state
                        IsotopicEnvelope putativeIsotopicEnvelope = FindIsotopicEnvelope(massIndex,
                            candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity,
                            testMostIntenseMass, chargeState, deconParams.DeconvolutionTolerancePpm, deconParams.IntensityRatioLimit,
                            monoisotopicMassPredictions);

                        if (putativeIsotopicEnvelope.Peaks.Count >= 2) //if there are at least two isotopes
                        {
                            //look for other charge states, using them for scoring and monoisotopic mass estimates
                            //need to use this method before comparing scores because it changes the score of the test envelope
                            int numOtherChargeStatesObserved = ObserveAdjacentChargeStates(putativeIsotopicEnvelope,
                                candidateForMostIntensePeakMz, massIndex, deconParams.DeconvolutionTolerancePpm,
                                deconParams.IntensityRatioLimit, deconParams.MinAssumedChargeState, deconParams.MaxAssumedChargeState,
                                monoisotopicMassPredictions);

                            //is this the best charge state for this peak?
                            if ((bestIsotopeEnvelopeForThisPeak == null ||
                                 putativeIsotopicEnvelope.Score >
                                 bestIsotopeEnvelopeForThisPeak
                                     .Score) && //and the score is better for this charge state than others
                                (putativeIsotopicEnvelope.Charge / 5 <=
                                 numOtherChargeStatesObserved)) //and if we suspect there to be multiple charge states and there are (higher the charge, more states expected, z=5, need 2 charge states, z=10, need 3 charge states, etc
                            {
                                putativeIsotopicEnvelope.SetMedianMonoisotopicMass(
                                    monoisotopicMassPredictions); //take the median mass from all of the isotopes (this is fine tuning!)
                                bestIsotopeEnvelopeForThisPeak = putativeIsotopicEnvelope;
                            }
                        }
                    }

                    if (bestIsotopeEnvelopeForThisPeak !=
                        null) //add this envelope (it might be wrong, but hopefully it has a low score and gets outscored later by the right thing)
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

        private IsotopicEnvelope FindIsotopicEnvelope(int massIndex, double candidateForMostIntensePeakMz, double candidateForMostIntensePeakIntensity, double testMostIntenseMass, int chargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, List<double> monoisotopicMassPredictions)
        {
            double[] theoreticalMasses = AverageResidueModel.GetAllTheoreticalMasses(massIndex);
            double[] theoreticalIntensities = AverageResidueModel.GetAllTheoreticalIntensities(massIndex);
            //add "most intense peak"
            var listOfObservedPeaks = new List<(double, double)> { (candidateForMostIntensePeakMz, candidateForMostIntensePeakIntensity) };
            var listOfRatios = new List<double> { theoreticalIntensities[0] / candidateForMostIntensePeakIntensity }; // theoreticalIntensities and theoreticalMasses are sorted by intensity, so first is most intense
            // Assuming the test peak is most intense...
            // Try to find the rest of the isotopes!
            double differenceBetweenTheorAndActualMass = testMostIntenseMass - theoreticalMasses[0]; //mass difference actual-theoretical for the tallest peak (not necessarily the monoisotopic)
            double totalIntensity = candidateForMostIntensePeakIntensity;
            double monoisotopicMass = testMostIntenseMass - AverageResidueModel.GetDiffToMonoisotopic(massIndex); //get the  monoisotopic by taking the most intense mass minus the expected mass difference between most intense and monoisotopic
            monoisotopicMassPredictions.Add(monoisotopicMass);
            for (int indexToLookAt = 1; indexToLookAt < theoreticalIntensities.Length; indexToLookAt++) //cycle through all theoretical peaks in this envelope from most intense to least intense
            {
                double theorMassThatTryingToFind = theoreticalMasses[indexToLookAt] + differenceBetweenTheorAndActualMass; //get the expected mass of the next most intense peak
                int matchedPeakIndex = spectrum.GetClosestPeakIndex(theorMassThatTryingToFind.ToMz(chargeState)); //find the experimental peak for that mass
                double matchedPeakMz = spectrum.XArray[matchedPeakIndex];
                double matchedPeakIntensity = spectrum.YArray[matchedPeakIndex];
                double matchedPeakMass = matchedPeakMz.ToMass(chargeState);

                //the nearest peak is out of tolerance, and it bounds every other peak's distance from the sought mass
                if (Math.Abs(matchedPeakMass - theorMassThatTryingToFind) / theorMassThatTryingToFind * 1e6 > deconvolutionTolerancePpm)
                {
                    break;
                }

                //the nearest peak is unusable, but issue #245: a usable isotope may sit just behind it, still in tolerance
                if (!Peak2satisfiesRatio(theoreticalIntensities[0], theoreticalIntensities[indexToLookAt], candidateForMostIntensePeakIntensity, matchedPeakIntensity, intensityRatioLimit)
                    || listOfObservedPeaks.Contains((matchedPeakMz, matchedPeakIntensity)))
                {
                    matchedPeakIndex = FindUsablePeakNearTheorMass(theorMassThatTryingToFind, chargeState, deconvolutionTolerancePpm,
                        theoreticalIntensities[0], theoreticalIntensities[indexToLookAt], candidateForMostIntensePeakIntensity,
                        intensityRatioLimit, listOfObservedPeaks);
                    if (matchedPeakIndex < 0)
                    {
                        break;
                    }

                    matchedPeakMz = spectrum.XArray[matchedPeakIndex];
                    matchedPeakIntensity = spectrum.YArray[matchedPeakIndex];
                    matchedPeakMass = matchedPeakMz.ToMass(chargeState);
                }

                //Found a match to an isotope peak for this charge state!
                listOfObservedPeaks.Add((matchedPeakMz, matchedPeakIntensity)); //add to observed list
                totalIntensity += matchedPeakIntensity; //add intensity
                listOfRatios.Add(theoreticalIntensities[indexToLookAt] / matchedPeakIntensity); //add ratio
                double monoisotopicMassFromThisPeak = monoisotopicMass + matchedPeakMass - theorMassThatTryingToFind;
                monoisotopicMassPredictions.Add(monoisotopicMassFromThisPeak);
            }

            return new IsotopicEnvelope(listOfObservedPeaks, monoisotopicMass, chargeState, totalIntensity, listOfRatios.StandardDeviation());
        }

        private int ObserveAdjacentChargeStates(IsotopicEnvelope originalEnvelope, double mostIntensePeakMz, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit, double minChargeToLookFor, double maxChargeToLookFor, List<double> monoisotopicMassPredictions)
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

        private bool FindChargeStateOfMass(IsotopicEnvelope originalEnvelope, int zToInvestigate, double mostAbundantNeutralIsotopeToInvestigate, int massIndex, double deconvolutionTolerancePpm, double intensityRatioLimit, List<double> monoisotopicMassPredictions)
        {
            //we know the mass and the charge that we're looking for, just see if the expected m/z and its isotopes are there or not
            double mostAbundantIsotopeMzForThisZTheoretical = mostAbundantNeutralIsotopeToInvestigate.ToMz(zToInvestigate);
            //let's go find that peak!
            int observedPeakIndex = spectrum.GetClosestPeakIndex(mostAbundantIsotopeMzForThisZTheoretical);

            double mostAbundantIsotopeMzObserved = spectrum.XArray[observedPeakIndex];
            double mostAbundantIsotopeMassObserved = mostAbundantIsotopeMzObserved.ToMass(zToInvestigate);
            //make sure the expected and observed peak are within the mass tolerance
            if (Math.Abs(mostAbundantIsotopeMassObserved - mostAbundantNeutralIsotopeToInvestigate) / mostAbundantNeutralIsotopeToInvestigate * 1e6 <= deconvolutionTolerancePpm)
            {
                //get the isotopic envelope for this peak and add the masses from all the peaks of the envelope to the monoisotopic mass predictions
                IsotopicEnvelope test = FindIsotopicEnvelope(massIndex, mostAbundantIsotopeMzObserved, spectrum.YArray[observedPeakIndex], mostAbundantIsotopeMassObserved,
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

        private (int start, int end) ExtractIndices(double minX, double maxX)
        {
            if (spectrum.XArray.Last() < minX || spectrum.XArray.First() > maxX)
                return (1, 0);

            return (
                spectrum.XArray.GetClosestIndex(minX, ArraySearchOption.Next),
                spectrum.XArray.GetClosestIndex(maxX, ArraySearchOption.Previous));
        }

        /// <summary>
        /// Nearest peak to <paramref name="theorMass"/> that is in tolerance, matches the predicted intensity and is
        /// unclaimed, or -1 if there is none. Called only when the closest peak itself is unusable (issue #245).
        /// </summary>
        private int FindUsablePeakNearTheorMass(double theorMass, int chargeState, double deconvolutionTolerancePpm,
            double peak1TheorIntensity, double peak2TheorIntensity, double peak1Intensity,
            double intensityRatioLimit, List<(double, double)> alreadyObservedPeaks)
        {
            int nearestIndex = spectrum.GetClosestPeakIndex(theorMass.ToMz(chargeState));

            //bound the window in the mass domain, since the acceptance test is a ppm check on neutral mass; converting
            //the other way would clip it by the proton offset
            double windowRadius = theorMass * deconvolutionTolerancePpm / 1e6;
            double minMz = (theorMass - windowRadius).ToMz(chargeState);
            double maxMz = (theorMass + windowRadius).ToMz(chargeState);

            int startIndex = nearestIndex;
            while (startIndex > 0 && spectrum.XArray[startIndex - 1] >= minMz)
            {
                startIndex--;
            }

            int endIndex = nearestIndex;
            while (endIndex < spectrum.XArray.Length - 1 && spectrum.XArray[endIndex + 1] <= maxMz)
            {
                endIndex++;
            }

            int bestIndex = -1;
            double bestMassError = double.PositiveInfinity;

            for (int index = startIndex; index <= endIndex; index++)
            {
                double candidateMz = spectrum.XArray[index];
                double candidateIntensity = spectrum.YArray[index];
                double massError = Math.Abs(candidateMz.ToMass(chargeState) - theorMass);

                //rank survivors by mass error: the ratio test has already screened out implausible abundances
                if (massError / theorMass * 1e6 > deconvolutionTolerancePpm
                    || candidateIntensity <= 0
                    || !Peak2satisfiesRatio(peak1TheorIntensity, peak2TheorIntensity, peak1Intensity, candidateIntensity, intensityRatioLimit)
                    || alreadyObservedPeaks.Contains((candidateMz, candidateIntensity)))
                {
                    continue;
                }

                if (massError < bestMassError)
                {
                    bestIndex = index;
                    bestMassError = massError;
                }
            }

            return bestIndex;
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
    }
}
