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

            //throw new NotImplementedException();

            // Because we've already generated all the isotopic envelopes we need, we should instead store their location by Protein
            // Need to link (PWSM, charge) keys to mz and intensity arrays 

            int numberOfBinsForIndexing = (int) ((SpectralParams.ScanRange.Maximum -
                                          SpectralParams.ScanRange.Minimum) *
                                          SpectralParams.BinsPerDalton).Ceiling(0);
            IndexedLibrarySpectra = new List<MinimalSpectrum>[numberOfBinsForIndexing];

            foreach (var keyValuePair in EnvelopeDictionary)
            {
                foreach (IsotopicEnvelope envelope in keyValuePair.Value)
                {
                    int binIndex =
                        (int)Math.Round(envelope.MostAbundantObservedIsotopicMass * SpectralParams.BinsPerDalton, 0);


                }
            }
        }

        private int GetBinIndex()
        {
            return 0;

        }
    }
}
