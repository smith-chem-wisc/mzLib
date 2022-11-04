using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Class for hosting deconvolution parameters common to all methods
    /// </summary>
    public abstract class DeconvolutionParameters
    {
        /// <summary>
        /// Constructor should initialize all fields that are used by every deconvolution algorithm
        /// </summary>
        public DeconvolutionParameters()
        {

        }
    }
}
