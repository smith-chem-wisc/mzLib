using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ExternalResults.BaseClasses
{
    public interface IQuantifiableResultFile
    {
        /// <summary>
        /// Returns every result in the result file as an IQuantifiableRecord
        /// </summary>
        /// <returns> Enumerable that contains identifications for a peptide </returns>
        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults();

        /// <summary>
        /// Links the file name associated with the protein to the raw file path of MassSpec data
        /// </summary>
        /// <param name="fileNames"> file name associated with each distinct record </param>
        /// <param name="filePath"> file path associated with each distinct record </param>
        /// <returns> Dictionary of file names and their associted full paths </returns>
        public Dictionary<string, string> FileNametoFilePath(List<string> fullFilePath);
    }
}
