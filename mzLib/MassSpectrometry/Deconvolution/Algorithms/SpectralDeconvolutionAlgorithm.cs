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
        public Dictionary<PeptideWithSetModifications, List<IsotopicEnvelope>> EnvelopeDictionary;
        private List<MinimalSpectrum>[] IndexedLibrarySpectra;
        // SpectrumIndexToPwsmMap maps the location of each spectrum within IndexedLibrarySpectra to its respective PeptideWithSetMods and charge
        private Dictionary<(int, int), (PeptideWithSetModifications, int charge)> SpectrumIndexToPwsmMap;
        public int MaxThreads; // This should maybe be in the Parameters abstract
        public SpectralDeconvolutionParameters SpectralParams;

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
                    IsotopicDistribution pwsmDistribution = IsotopicDistribution.GetDistribution(pwsm.FullChemicalFormula);
                    for (int charge = SpectralParams.MinAssumedChargeState;
                         charge <= SpectralParams.MinAssumedChargeState;
                         charge++)
                    {
                        double theoreticalMz = pwsm.MonoisotopicMass.ToMz(charge);
                        if (SpectralParams.ScanRange.Contains(theoreticalMz))
                        {
                            EnvelopeDictionary[pwsm].Add(new IsotopicEnvelope(pwsmDistribution, charge));
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

            int numberOfBinsForIndexing = (int) ((SpectralParams.ScanRange.Maximum -
                                          SpectralParams.ScanRange.Minimum) *
                                          SpectralParams.BinsPerDalton).Ceiling(0);
            IndexedLibrarySpectra = new List<MinimalSpectrum>[numberOfBinsForIndexing];
            SpectrumIndexToPwsmMap = new();

            foreach (var keyValuePair in EnvelopeDictionary)
            {
                foreach (IsotopicEnvelope envelope in keyValuePair.Value)
                {
                    int binIndex = (int)Math.Round(envelope.MostAbundantObservedIsotopicMass * SpectralParams.BinsPerDalton, 0);

                    MinimalSpectrum envelopeMinimalSpectrum = new MinimalSpectrum(envelope.MzArray, envelope.IntensityArray);
                    IndexedLibrarySpectra[binIndex].Add(envelopeMinimalSpectrum);
                    SpectrumIndexToPwsmMap.Add(
                        (binIndex, IndexedLibrarySpectra[binIndex].Count - 1), // tuple consisting of bin index and list position of MinimalSpectrum object
                        (keyValuePair.Key, envelope.Charge) // tuple consisting of PeptideWithSetMods and charge state
                        );

                    // In situations where the most abundant isotope frequency is close to the second most abundant isotope's frequency
                    // ( ratio >= IsotopicEnvelope.AmbiguityRatioMinimum),
                    // The Spectrum is stored in the index of the second most abundant isotope as well
                    if(envelope.SecondMostAbundantObservedIsotopicMass > 0)
                    {
                        binIndex = (int)Math.Round((double)envelope.SecondMostAbundantObservedIsotopicMass * SpectralParams.BinsPerDalton, 0);
                        IndexedLibrarySpectra[binIndex].Add(envelopeMinimalSpectrum);
                        SpectrumIndexToPwsmMap.Add(
                            (binIndex, IndexedLibrarySpectra[binIndex].Count - 1),
                            (keyValuePair.Key, envelope.Charge) 
                            );
                    }
                }
            }
        }

    }
}
