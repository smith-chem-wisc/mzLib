using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    /// <summary>
    /// Outlines behavior to turn results into an IEnumerable of IQuantifiableRecords 
    /// and to create the dictionary linking file names from the external result files 
    /// to their local file paths which are used to make the identification object
    /// </summary>
    public interface IQuantifiableResultFile : IResultFile
    {
        /// <summary>
        /// Returns every result in the result file as an IQuantifiableRecord
        /// </summary>
        /// <returns> Enumerable that contains identifications for a peptide </returns>
        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults();

        /// <summary>
        /// Links the file name associated with the protein to the raw file path of MassSpec data
        /// </summary>
        /// <param name="fullFilePath"> list of file paths associated with each distinct record </param>
        /// <returns> Dictionary of file names and their associted full paths </returns>
        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath);
    }
}