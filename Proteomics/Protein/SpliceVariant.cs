using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class SpliceVariant
    {
        public SpliceVariant(string accession, string name, string type, List<SequenceVariation> spliceVariations)
        {
            Accession = accession;
            Name = name;
            Type = type;
            SpliceVariations = spliceVariations;
        }

        public string Accession { get; private set; }
        public string Name { get; private set; }
        public string Type { get; private set; } // "displayed", "described", "not described"
        public List<SequenceVariation> SpliceVariations { get; private set; }

        public string GetDescription()
        {
            // missing 1 to 8
            throw new NotImplementedException();
        }
    }
}
