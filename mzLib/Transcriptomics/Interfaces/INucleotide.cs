using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Transcriptomics.Modifications;

namespace Transcriptomics
{
    public interface INucleotide : IHasChemicalFormula
    {
        char Letter { get; }
        string Symbol { get; }
        ModificationSite ModificationSite { get; }
    }
}
