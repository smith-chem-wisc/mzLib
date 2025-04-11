using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Parent class defining minimum requirement to be used <see cref="Deconvoluter"/> 
    /// </summary>
    public abstract class DeconvolutionAlgorithm
    {
        // For ClassicDeconv. If not used elsewhere, move to that class
        #region Averagine Stuff

        protected const int numAveraginesToGenerate = 1500;
        protected static readonly double[][] allMasses = new double[numAveraginesToGenerate][];
        protected static readonly double[][] allIntensities = new double[numAveraginesToGenerate][];
        protected static readonly double[] mostIntenseMasses = new double[numAveraginesToGenerate];
        protected static readonly double[] diffToMonoisotopic = new double[numAveraginesToGenerate];

        static DeconvolutionAlgorithm()
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

        #endregion

        protected readonly DeconvolutionParameters DeconvolutionParameters;

        /// <summary>
        /// Constructor for deconvolution algorithms, nothing should be added to child constructors
        /// </summary>
        /// <param name="deconParameters">parameters to use for deconvolution</param>
        /// <exception cref="MzLibException">thrown in deconvolution parameters did not instantiate fields required by algorithm</exception>
        protected DeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
        {
            DeconvolutionParameters = deconParameters;
        }

        /// <summary>
        /// Deconvolutes a mass spectrum
        /// </summary>
        /// <param name="spectrum">spectrum to be deconvoluted</param>
        /// <param name="range">Range of peaks to deconvolute</param>
        /// <returns></returns>
        internal abstract IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range);
    }
}
