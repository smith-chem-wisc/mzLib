using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    public class QuantifiedProteinGroup
    {
        public string Name { get; set; }
        public Dictionary<string, QuantifiedProtein> Proteins { get; set; }

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
