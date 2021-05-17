using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    /// <summary>
    /// Store information pertaining each unique detected protein
    /// </summary>
    public class ProteinRowInfo
    {
        public ProteinRowInfo()
        {
            SamplesIntensityData = new Dictionary<string, double>();
        }

        /// <summary>
        /// Stores detected intensity values of protein in each sample. 
        /// Maps condition's sample name to the intensity value
        /// </summary>
        public Dictionary<string, double> SamplesIntensityData { get; set; }

        /// <summary>
        /// Unique protein ID
        /// </summary>
        public string ProteinID { get; set; }
    }
}