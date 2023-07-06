using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace Transcriptomics.Fragmentation
{
    public class ChemicalFormulaFragment : Fragment, IHasChemicalFormula
    {
        public ChemicalFormulaFragment(FragmentType type, int number, string chemicalFormula, NucleicAcid parent) :
            this(type, number, ChemicalFormula.ParseFormula(chemicalFormula), parent)
        {
        }

        public ChemicalFormulaFragment(FragmentType type, int number, ChemicalFormula chemicalFormula, NucleicAcid parent)
            : base(type, number, chemicalFormula.MonoisotopicMass, parent)
        {
            ThisChemicalFormula = chemicalFormula;
        }

        public ChemicalFormula ThisChemicalFormula { get; }
    }
}
