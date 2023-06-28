using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace Transcriptomics
{
    public class RNA : NucleicAcid
    {
        public RNA(string sequence) : base(sequence)
        {
        }

        public RNA(string sequence, IHasChemicalFormula fivePrimeTerm, IHasChemicalFormula threePrimeTerm) : base(sequence, fivePrimeTerm, threePrimeTerm)
        {
        }
    }
}
