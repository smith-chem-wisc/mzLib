using System.Collections.Generic;
using Proteomics;

namespace UsefulProteomicsDatabases
{
    public class Modification
    {
        private IEnumerable<string> motifs;
        private string uniprotAC;
        private Dictionary<string, HashSet<string>> uniprotDR;
        private string uniprotID;
        private ModificationSites uniprotPP;
        private ModificationSites uniprotTG;

        public Modification(string uniprotID, string uniprotAC, ModificationSites uniprotTG, ModificationSites uniprotPP, Dictionary<string, HashSet<string>> uniprotDR, IEnumerable<string> motifs)
        {
            this.uniprotID = uniprotID;
            this.uniprotAC = uniprotAC;
            this.uniprotTG = uniprotTG;
            this.uniprotPP = uniprotPP;
            this.uniprotDR = uniprotDR;
            this.motifs = motifs;
        }
    }
}