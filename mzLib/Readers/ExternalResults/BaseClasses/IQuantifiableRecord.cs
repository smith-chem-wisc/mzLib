using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ExternalResults.BaseClasses
{
    public interface IQuantifiableRecord
    {
        /// <summary>
        /// The full file path to the MS Data file in which the identification was made
        /// </summary>
        public string FullFilePath { get; set; }

        /// <summary>
        /// A list of tuples, each of which represent a protein. Each tuple contains the accession number, gene name, and organism associated with the given result
        /// </summary>
        public List<(string proteinAccessions, string geneName, string organism)> proteinGroupInfos { get; set; }

        /// <summary>
        /// The amino acid sequence of the identified peptide
        /// </summary>
        public string BaseSequence { get; set; }

        /// <summary>
        /// The amino acid sequence and the associated post-translation modifications of the identified peptide
        /// </summary>
        public string ModifiedSequence { get; set; }

        /// <summary>
        /// the retention time (in minutes) associated with the result
        /// </summary>
        public double RetentionTime { get; set; }

        /// <summary>
        /// The charge state associated with the result
        /// </summary>
        public int ChargeState { get; set; }

        /// <summary>
        /// Defines whether or not the result is a decoy identification
        /// </summary>
        public bool IsDecoy { get; set; }

        /// <summary>
        /// The mass of the monoisotopic peptide (i.e., no c13 or n15 atoms are present, the lowest possible mass)
        /// </summary>
        public double MonoisotopicMass { get; set; }

     
    }
}
