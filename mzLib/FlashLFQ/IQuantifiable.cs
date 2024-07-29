using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ 
{ 
    internal interface IQuantifiable
    {
        /// <summary>
        /// A dictionary that links the name of a data file to the associated SpectraFileInfo object
        /// </summary>
        public Dictionary<string, SpectraFileInfo> SpectraFileDict { get; set; }

        /// <summary>
        /// A dictionary that links protein accessions to ProteinGroup objects
        /// </summary>
        public Dictionary<string, ProteinGroup> ProteinGroupDict { get; set; }

        /// <summary>
        /// Returns each result in the ResultFile as an Identification object
        /// </summary>
        /// <param name="dataFileDirectory"> a local file path to the directory containing all the ms data files that were searched </param>
        /// <returns> enumerable that contains identifications for a peptide </returns>
        public IEnumerable<Identification> GetIdentifications(string dataFileDirectory); 

    }
}
