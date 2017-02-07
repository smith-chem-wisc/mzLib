using System.Collections.Generic;
using Proteomics;
using UsefulProteomicsDatabases.Generated;
using System;

namespace UsefulProteomicsDatabases
{
    public class Modification
    {

        private readonly Tuple<string, string> ac;
        private readonly Dictionary<string, HashSet<string>> linksToOtherDbs;
        private readonly string id;
        private readonly ModificationSites position;
        private readonly string site;

        private static readonly Dictionary<position_t, ModificationSites> positionDict = new Dictionary<position_t, ModificationSites>
            {
            {position_t.AnyCterm, ModificationSites.PepC },
            {position_t.ProteinCterm, ModificationSites.ProtC },
            {position_t.Anywhere, ModificationSites.Any },
            {position_t.AnyNterm, ModificationSites.NPep },
            {position_t.ProteinNterm, ModificationSites.NProt }
            };

        public Modification(string id, Tuple<string, string> unimodAC, string tg, position_t pos)
        {
            this.id = id;
            this.ac = unimodAC;
            this.site = tg;
            this.position = positionDict[pos];
        }

        public Modification(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, Dictionary<string, HashSet<string>> uniprotDR)
        {
            this.id = uniprotID;
            this.ac = uniprotAC;
            this.site = uniprotTG;
            this.position = uniprotPP;
            this.linksToOtherDbs = uniprotDR;
        }
    }
}