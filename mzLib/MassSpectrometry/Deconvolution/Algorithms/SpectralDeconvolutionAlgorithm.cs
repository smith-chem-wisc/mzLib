using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry.Deconvolution;
using MassSpectrometry.Deconvolution.Scoring;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    public class SpectralDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        // TODO: Make a charge state envelope class, complete with "MostAbundantChargeState"
        public Dictionary<PeptideWithSetModifications, List<IsotopicEnvelope>> EnvelopeDictionary { get; private set; }

        // Consider defining this as a jagged array to increase performance
        public List<MinimalSpectrum>[,] IndexedLibrarySpectra { get; private set; }
        // SpectrumIndexToPwsmMap maps the location of each spectrum within IndexedLibrarySpectra to its respective PeptideWithSetMods and charge
        public Dictionary<(int, int, int), PeptideWithSetModifications> SpectrumIndexToPwsmMap { get; private set; }
        public int MaxThreads; // This should be in the Parameters abstract 
        public SpectralDeconvolutionParameters SpectralParams { get; }
        public PpmTolerance PpmTolerance { get; }
        public Scorer Scorer { get; }

        public SpectralDeconvolutionAlgorithm(DeconvolutionParameters parameters) : base(parameters)
        {
            var deconvolutionParameters = DeconvolutionParameters as SpectralDeconvolutionParameters;
            if (deconvolutionParameters == null)
            {
                throw new MzLibException(
                    "Improper Deconvolution Parameters were pass to the SpectralDeconvolutionAlgorithm");
            }
            else
            {
                SpectralParams = deconvolutionParameters;
            }

            PpmTolerance = new PpmTolerance(parameters.DeconvolutionTolerancePpm);

            FindLibraryEnvelopes();
            IndexEnvelopes();
            Scorer = new Scorer(Scorer.ScoringMethods.SpectralContrastAngle, PpmTolerance);

        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum)
        {

            if (spectrum == null || spectrum.Size == 0)
            {
                yield break;
            }

            // For each charge state (key) store the indices corresponding to every potential isotopic envelope (value)
            Dictionary<int, List<MinimalSpectrum>> potentialEnvelopes = FindPotentialEnvelopes(spectrum);

            // iterate through charge states (potentially not necessary/performant. Could flatten)
            foreach (var keyValuePair in potentialEnvelopes)
            {
                int chargeBinIndex = keyValuePair.Key - SpectralParams.MinAssumedChargeState;
                // iterate through potential envelopes
                foreach (var experimentalSpectrum in keyValuePair.Value)
                {
                    double mostAbundantMz =  experimentalSpectrum.MostAbundantMz;
                    int massBinIndex = (int)Math.Floor((mostAbundantMz - SpectralParams.ScanRange.Minimum) *
                                                   SpectralParams.BinsPerDalton);
                    if (!IndexedLibrarySpectra[massBinIndex, chargeBinIndex].IsNotNullOrEmpty()) continue; // continue if there are no corresponding library spectra

                    int? bestMatchListPosition = null;
                    int currentListPosition = 0;
                    double bestFoundScore = Scorer.PoorScore;
                    // Score against matching theoretical envelopes
                    foreach (MinimalSpectrum theoreticalSpectrum in IndexedLibrarySpectra[massBinIndex, chargeBinIndex])
                    {
                        // TODO: Rename to FindBestScore
                        if (Scorer.TestForScoreImprovement(
                                Scorer.Score(experimentalSpectrum,theoreticalSpectrum),
                                bestFoundScore,
                                out double betterScore)
                            )
                        {
                            bestMatchListPosition = currentListPosition;
                            bestFoundScore = betterScore;
                        }
                        currentListPosition++;
                    }

                    if (bestMatchListPosition.HasValue && 
                        SpectrumIndexToPwsmMap.TryGetValue((massBinIndex, chargeBinIndex, (int)bestMatchListPosition), out var pwsmMatch))
                    {
                        yield return new IsotopicEnvelope(experimentalSpectrum, pwsmMatch, bestFoundScore);
                    }
                    else
                    {
                        //TODO: Add some averagine bullshit here
                    }
                    
                }
            }
        }

        /// <summary>
        /// Populates the EnvelopeDictionary by digesting each protein in the parameters into PeptideWithSetMods,
        /// then calculating an isotopic envelope for each charge state from min to max assumed charge state
        /// </summary>
        private void FindLibraryEnvelopes()
        {
            EnvelopeDictionary = new();
            

            //TODO: Parallelize this section of the code
            foreach (Protein protein in SpectralParams.Proteins)
            {
                // I'm not sure if calling protein.Digest within the foreach statement would call the method anew for every loop
                IEnumerable<PeptideWithSetModifications> uniquePeptides = protein.Digest(
                    SpectralParams.DigestionParams, SpectralParams.FixedModifications,
                    SpectralParams.VariableModifications, SpectralParams.SilacLabels,
                    topDownTruncationSearch: SpectralParams.FindTopDownTruncationProducts);

                foreach (PeptideWithSetModifications pwsm in uniquePeptides)
                {
                    EnvelopeDictionary.Add(pwsm, new List<IsotopicEnvelope>());
                    IsotopicDistribution pwsmDistribution = IsotopicDistribution.GetDistribution(pwsm.FullChemicalFormula,
                        fineResolution: SpectralParams.FineResolutionForIsotopicDistribution,
                        minProbability: SpectralParams.MinProbabilityForIsotopicDistribution);

                    // iterates through all possible charge states, from largest to smallest.
                    // Any isotopic envelope whose most abundant peak would fall within the scan range is written to the envelope dictionary
                    // Once the mass to charge ratio of the most abundant peak is greater than the scan range maximum, the loop breaks
                    for (int charge = SpectralParams.MaxAssumedChargeState;
                         charge >= SpectralParams.MinAssumedChargeState;
                         charge--)
                    {
                        double theoreticalMz = pwsm.MostAbundantMass.ToMz(charge);
                        if (SpectralParams.ScanRange.Contains(theoreticalMz))
                        {
                            EnvelopeDictionary[pwsm].Add(new 
                                IsotopicEnvelope(pwsmDistribution, charge, SpectralParams.AmbiguityThresholdForIsotopicDistribution));
                        }
                        else if (SpectralParams.ScanRange.CompareTo(theoreticalMz) < 0)
                        {
                            break;
                        }
                    }
                }

            }
        }

        /// <summary>
        /// For each envelope in Envelope Dictionary, indexes it according to mass and charge,
        /// resulting in a 2D array of lists of minimal spectra
        /// </summary>
        private void IndexEnvelopes()
        {

            int numberOfBinsForIndexing = (int) (SpectralParams.ScanRange.Width * SpectralParams.BinsPerDalton).Ceiling(0);
            IndexedLibrarySpectra = new List<MinimalSpectrum>[numberOfBinsForIndexing,
                SpectralParams.MaxAssumedChargeState + 1 - SpectralParams.MinAssumedChargeState];
            SpectrumIndexToPwsmMap = new();

            foreach (var keyValuePair in EnvelopeDictionary)
            {
                foreach (IsotopicEnvelope envelope in keyValuePair.Value)
                {
                    int massBinIndex = (int)Math.Floor((envelope.MostAbundantObservedIsotopicMz - SpectralParams.ScanRange.Minimum) * 
                                                     SpectralParams.BinsPerDalton);
                    int chargeBinIndex = envelope.Charge - SpectralParams.MinAssumedChargeState;
                    if (IndexedLibrarySpectra[massBinIndex, chargeBinIndex] == null)
                    {
                        IndexedLibrarySpectra[massBinIndex, chargeBinIndex] = new();
                    }
                    MinimalSpectrum envelopeMinimalSpectrum = new MinimalSpectrum(envelope.MzArray, envelope.IntensityArray, envelope.Charge);
                    IndexedLibrarySpectra[massBinIndex, chargeBinIndex].Add(envelopeMinimalSpectrum);
                    SpectrumIndexToPwsmMap.Add(
                        (massBinIndex, chargeBinIndex, IndexedLibrarySpectra[massBinIndex, chargeBinIndex].Count - 1), // tuple consisting of bin index (mass, charge) and list position of MinimalSpectrum object
                        keyValuePair.Key // tuple consisting of PeptideWithSetMods and charge state
                        );

                    // In situations where the most abundant isotope frequency is close to the second most abundant isotope's frequency
                    // ( ratio >= IsotopicEnvelope.AmbiguityRatioMinimum), 
                    // The Spectrum is stored in the index of the second most abundant isotope as well
                    if(envelope.SecondMostAbundantObservedIsotopicMz > 0 )
                    {
                        // Ceiling or floor????
                        int secondBinIndex = (int)Math.Floor( 
                            ((double)envelope.SecondMostAbundantObservedIsotopicMz - SpectralParams.ScanRange.Minimum ) * SpectralParams.BinsPerDalton);
                        if (secondBinIndex != massBinIndex)
                        {
                            if (IndexedLibrarySpectra[secondBinIndex, chargeBinIndex] == null) IndexedLibrarySpectra[secondBinIndex, chargeBinIndex] = new();
                            IndexedLibrarySpectra[secondBinIndex, chargeBinIndex].Add(envelopeMinimalSpectrum);
                            SpectrumIndexToPwsmMap.Add(
                                (secondBinIndex, chargeBinIndex, IndexedLibrarySpectra[secondBinIndex, chargeBinIndex].Count - 1),
                                keyValuePair.Key
                            );
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Iterates through all peaks in a spectrum to find all potential isotopic envelopes.
        /// It does this by examining the spacing of peaks in the m/z domain
        /// e.g. for charge of 2, a peak at 200 m/z would result in a search for a peak at 200.5 and 201 m/z
        ///      if either is found, the process continues until SpectralParams.MaxConsecutiveMissedIsotopicPeaks number of consecutive
        ///      isotope peaks are missed
        /// Anything consistent with an isotopic envelope in a given charge state is stored in the dictionary
        /// </summary>
        /// <param name="spectrum"></param>
        /// <returns></returns>
        private Dictionary<int, List<MinimalSpectrum>> FindPotentialEnvelopes(MzSpectrum spectrum)
        {

            // For each charge state (key) store the indices corresponding to every potential isotopic envelope (value)
            Dictionary<int, List<MinimalSpectrum>> potentialEnvelopes = new();

            for (int charge = SpectralParams.MinAssumedChargeState; charge <= SpectralParams.MaxAssumedChargeState; charge++)
            {
                List<int> indicesOfKnownPeaks = new();

                // Spectrum Search Loop
                for (int i = 0; i < spectrum.Size; i++)
                {
                    if (indicesOfKnownPeaks.Contains(i))
                    {
                        continue;
                    }
                    List<int> envelopeIndices = new();
                    envelopeIndices.Add(i);

                    // Envelope Search Loop
                    for (int j = i + 1; j < spectrum.Size; j++)
                    {
                        if (PpmTolerance.Within(spectrum.XArray[j],
                                spectrum.XArray[envelopeIndices.Last()] + Constants.C13MinusC12 / charge))
                        {
                            envelopeIndices.Add(j);
                        }
                        else if (spectrum.XArray[j] > PpmTolerance.GetMaximumValue(spectrum.XArray[envelopeIndices.Last()] +
                                     (1 + SpectralParams.MaxConsecutiveMissedIsotopicPeaks) * Constants.C13MinusC12 / charge))
                        {
                            // exit the Envelope loop if we missed more consecutive isotopic peaks than were allowed
                            break;
                        }
                    }

                    // Convert to MinimalSpectrum here? Write helper function to do so?
                    if (envelopeIndices.Count > 1)
                    {
                        if (!potentialEnvelopes.ContainsKey(charge)) potentialEnvelopes.Add(charge, new());
                        potentialEnvelopes[charge].Add(GetMinimalSpectrumFromIndices(spectrum, envelopeIndices, charge));
                        indicesOfKnownPeaks.AddRange(envelopeIndices);
                    }
                }
            }

            return potentialEnvelopes;
        }


        private static MinimalSpectrum GetMinimalSpectrumFromIndices(MzSpectrum spectrum, List<int> indices, int charge = 0)
        {
            double[] mzArray = new double[indices.Count];
            double[] intensityArray = new double[indices.Count];
            for (int i = 0; i < indices.Count; i++)
            {
                mzArray[i] = spectrum.XArray[indices[i]];
                intensityArray[i] = spectrum.YArray[indices[i]];
            }

            return new MinimalSpectrum(mzArray, intensityArray, charge);
        }

    }
}
