using Chemistry;
using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public enum LabelType { Variable, Fixed }

    public class ChemicalLabel
    {
        public string Label;
        public ChemicalFormula LabelChemicalFormula;
        public double LabelMass;
        public List<char> LabelResidues;
        public LabelType LabelType;

        public ChemicalLabel(ChemicalFormula formula)
        {

        }

        public static Dictionary<string, List<ChemicalLabel>> KnownChemicalLabels = new Dictionary<string, List<ChemicalLabel>>
        {
            { "K(+8.048)", new List<ChemicalLabel> { } },
            { "TMT-11", new List<ChemicalLabel> 
            { 
                new ChemicalLabel(ChemicalFormula.ParseFormula("C8N1H15")) } 
            },
        };
    }
}
