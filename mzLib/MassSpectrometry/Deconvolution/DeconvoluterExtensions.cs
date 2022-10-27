using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Class to provide extensions for other objects besides MsDataScan to be deconvoluted
    /// Methods can be generated for each deconvolution type, passing parameters that would otherwise be found in the MsDataScan
    /// </summary>
    public static class DeconvoluterExtensions
    {
        public static IEnumerable<IsotopicEnvelope> ClassicDeconvoluteMzSpectra(this Deconvoluter deconvoluter,
            MzSpectrum spectrum, MzRange range)
        {
            if (deconvoluter.DeconvolutionType != DeconvolutionTypes.ClassicDeconvolution)
            {
                throw new MzLibException("Deconvoluter is not of correct type for this extension method");
            }
            else
            {
                ((ClassicDeconvolutionParameters)deconvoluter.DeconvolutionParameters).Range = range;
                return deconvoluter.DeconvolutionAlgorithm.Deconvolute(spectrum);
            }
        }
    }
}
