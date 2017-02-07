using System.Collections.Generic;
using Proteomics;
using UsefulProteomicsDatabases.Generated;
using System;

namespace UsefulProteomicsDatabases
{
    internal class ModificationWithMass : Modification
    {
        private readonly double neutralLoss;
        private readonly IEnumerable<double> massesObserved;
        private readonly double monoisotopicMass;
        private readonly IEnumerable<double> diagnosticIons;


        public ModificationWithMass(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotDR)
        {
            this.monoisotopicMass = uniprotMM;
            this.neutralLoss = neutralLoss;
            this.massesObserved = massesObserved;
            this.diagnosticIons = diagnosticIons;
        }

        public ModificationWithMass(string id, Tuple<string,string> ac, string tg, position_t pos, double mm, double neutralLoss)
            : base(id, ac, tg, pos)
        {
            this.neutralLoss = neutralLoss;
            this.monoisotopicMass = mm;
            massesObserved = new HashSet<double> { mm };
        }
    }
}