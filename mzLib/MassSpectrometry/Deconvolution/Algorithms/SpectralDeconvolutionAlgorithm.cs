using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry.Deconvolution.Parameters;
using MassSpectrometry.Deconvolution;
using MassSpectrometry.Proteomics;
using MassSpectrometry.Proteomics.ProteolyticDigestion;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    public class SpectralDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        // TODO: Make a charge state envelope class, complete with "MostAbundantChargeState"
        public Dictionary<PeptideWithSetModifications, List<IsotopicEnvelope>> EnvelopeDictionary { get; private set; }

        // Consider defining this as a jagged array to increase performance
        public List<MinimalSpectrum>[] IndexedLibrarySpectra { get; private set; }
        // SpectrumIndexToPwsmMap maps the location of each spectrum within IndexedLibrarySpectra to its respective PeptideWithSetMods and charge
        public Dictionary<(int, int), (PeptideWithSetModifications pwsm, int charge)> SpectrumIndexToPwsmMap { get; private set; }
        public int MaxThreads; // This should be in the Parameters abstract {
        public SpectralDeconvolutionParameters SpectralParams { get; }

        public SpectralDeconvolutionAlgorithm(SpectralDeconvolutionParameters parameters) : base(parameters)
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

            // Calculate every species you expect to see
            FindLibraryEnvelopes();
            // Index envelopes
            IndexEnvelopes();
        }

        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum)
        {
            throw new NotImplementedException();
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

        private void IndexEnvelopes()
        {

            int numberOfBinsForIndexing = (int) (SpectralParams.ScanRange.Width * SpectralParams.BinsPerDalton).Ceiling(0);
            IndexedLibrarySpectra = new List<MinimalSpectrum>[numberOfBinsForIndexing];
            SpectrumIndexToPwsmMap = new();

            foreach (var keyValuePair in EnvelopeDictionary)
            {
                foreach (IsotopicEnvelope envelope in keyValuePair.Value)
                {
                    int binIndex = (int)Math.Floor((envelope.MostAbundantObservedIsotopicMz - SpectralParams.ScanRange.Minimum) * 
                                                     SpectralParams.BinsPerDalton);
                    if (IndexedLibrarySpectra[binIndex] == null) IndexedLibrarySpectra[binIndex] = new();
                    MinimalSpectrum envelopeMinimalSpectrum = new MinimalSpectrum(envelope.MzArray, envelope.IntensityArray);
                    IndexedLibrarySpectra[binIndex].Add(envelopeMinimalSpectrum);
                    SpectrumIndexToPwsmMap.Add(
                        (binIndex, IndexedLibrarySpectra[binIndex].Count - 1), // tuple consisting of bin index and list position of MinimalSpectrum object
                        (keyValuePair.Key, envelope.Charge) // tuple consisting of PeptideWithSetMods and charge state
                        );

                    // In situations where the most abundant isotope frequency is close to the second most abundant isotope's frequency
                    // ( ratio >= IsotopicEnvelope.AmbiguityRatioMinimum), 
                    // The Spectrum is stored in the index of the second most abundant isotope as well
                    if(envelope.SecondMostAbundantObservedIsotopicMz > 0 )
                    {
                        // Ceiling or floor????
                        int secondBinIndex = (int)Math.Floor( 
                            ((double)envelope.SecondMostAbundantObservedIsotopicMz - SpectralParams.ScanRange.Minimum ) * SpectralParams.BinsPerDalton);
                        if (secondBinIndex != binIndex)
                        {
                            if (IndexedLibrarySpectra[secondBinIndex] == null) IndexedLibrarySpectra[secondBinIndex] = new();
                            IndexedLibrarySpectra[secondBinIndex].Add(envelopeMinimalSpectrum);
                            SpectrumIndexToPwsmMap.Add(
                                (secondBinIndex, IndexedLibrarySpectra[secondBinIndex].Count - 1),
                                (keyValuePair.Key, envelope.Charge)
                            );
                        }
                    }
                }
            }
        }

    }
}
