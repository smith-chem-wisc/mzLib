using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry.Deconvolution.Parameters;
using MassSpectrometry.Proteomics;
using MassSpectrometry.Proteomics.ProteolyticDigestion;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    public class SpectralDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        // TODO: Make a charge state envelope class, complete with "MostAbundantChargeState"
        public Dictionary<PeptideWithSetModifications, List<IsotopicEnvelope>> EnvelopeDictionary;
        public List<Protein> ProteinsInLibrary;
        public int MaxThreads;

        public SpectralDeconvolutionAlgorithm(SpectralDeconvolutionParameters parameters) : base(parameters)
        {
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
            var deconvolutionParameters = DeconvolutionParameters as SpectralDeconvolutionParameters;
            if (deconvolutionParameters == null)
            {
                throw new MzLibException(
                    "Improper Deconvolution Parameters were pass to the SpectralDeconvolutionAlgorithm");
            } 

            //TODO: Parallelize this section of the code
            foreach (Protein protein in deconvolutionParameters.Proteins)
            {
                // I'm not sure if calling protein.Digest within the foreach statement would call the method anew for every loop
                IEnumerable<PeptideWithSetModifications> uniquePeptides = protein.Digest(
                    deconvolutionParameters.DigestionParams, deconvolutionParameters.FixedModifications,
                    deconvolutionParameters.VariableModifications, deconvolutionParameters.SilacLabels,
                    topDownTruncationSearch: deconvolutionParameters.FindTopDownTruncationProducts);

                foreach (PeptideWithSetModifications pwsm in uniquePeptides)
                {
                    EnvelopeDictionary.Add(pwsm, new List<IsotopicEnvelope>());
                    IsotopicDistribution pwsmDistribution = IsotopicDistribution.GetDistribution(pwsm.FullChemicalFormula);
                    for (int charge = deconvolutionParameters.MinAssumedChargeState;
                         charge <= deconvolutionParameters.MinAssumedChargeState;
                         charge++)
                    {
                        double theoreticalMz = pwsm.MonoisotopicMass.ToMz(charge);
                        if (deconvolutionParameters.ScanRange.Contains(theoreticalMz))
                        {
                            EnvelopeDictionary[pwsm].Add(new IsotopicEnvelope(pwsmDistribution, charge));
                        }
                        else if (deconvolutionParameters.ScanRange.CompareTo(theoreticalMz) < 0)
                        {
                            break;
                        }
                    }
                }

            }
        }

        private void IndexEnvelopes()
        {
            throw new NotImplementedException();
            // Because we've already generated all the isotopic envelopes we need, we should instead store their location by (Protein
            // Need to link (PWSM, charge) keys to mz and intensity arrays 
        }
    }
}
