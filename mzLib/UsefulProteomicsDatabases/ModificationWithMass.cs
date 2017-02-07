using System.Collections.Generic;
using Proteomics;

namespace UsefulProteomicsDatabases
{
    internal class ModificationWithMass : Modification
    {
        private IEnumerable<double> massesOnFragments;
        private IEnumerable<double> massesObserved;
        private double? uniprotMM;

        public ModificationWithMass(string uniprotID, string uniprotAC, ModificationSites uniprotTG, ModificationSites uniprotPP, double? uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, IEnumerable<double> massesOnFragments, IEnumerable<double> massesObserved, IEnumerable<string> motifs)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotDR, motifs)
        {
            this.uniprotMM = uniprotMM;
            this.massesOnFragments = massesOnFragments;
            this.massesObserved = massesObserved;
        }
    }
}