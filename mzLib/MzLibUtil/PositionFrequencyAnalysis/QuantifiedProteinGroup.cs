using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// Represents a group of proteins for quantification purposes.
    /// </summary>
    public class QuantifiedProteinGroup
    {
        /// <summary>
        /// The name of the protein group, typically a concatenation of protein accessions in the 
        /// format "ProteinA;ProteinB", "ProteinA|ProteinB", or "ProteinA;ProteinB|ProteinC".
        /// </summary>
        public string Name { get; set; }

        public string GeneName { get; set; }
        public string Organism { get; set; }

        /// <summary>
        /// Dictionary mapping protein accessions to their corresponding QuantifiedProtein objects.
        /// </summary>

        public Dictionary<string, QuantifiedProtein> Proteins { get; set; }

        /// <summary>
        /// Initializes a new protein group with the specified name and optional proteins.
        /// </summary>
        public QuantifiedProteinGroup(string name, Dictionary<string, QuantifiedProtein> proteins = null)
        {
            proteins = proteins ?? new Dictionary<string, QuantifiedProtein>();
            string splitPattern = @";|\|";
            var proteinAccessions = Regex.Split(name, splitPattern);
            if (proteinAccessions.Length == proteins.Count && proteinAccessions.OrderBy(x => x).SequenceEqual(proteins.Keys.OrderBy(x => x)) || proteins.IsNullOrEmpty())
            {
                Name = name;
                Proteins = proteins ?? new Dictionary<string, QuantifiedProtein>();
            }
            else
            {
                throw new Exception("The number of proteins provided does not match the number of proteins in the protein group name.");
            }
        }
    }
}
