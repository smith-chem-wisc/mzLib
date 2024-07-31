using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ExternalResults.BaseClasses
{
    public interface IQuantifiable
    {
        /// <summary>
        /// Returns every result in the result file as an IQuantifiableRecord
        /// </summary>
        /// <returns> Enumerable that contains identifications for a peptide </returns>
        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults();

        public Dictionary<string, string> FileNametoFilePath { get; set; }
    }
}
