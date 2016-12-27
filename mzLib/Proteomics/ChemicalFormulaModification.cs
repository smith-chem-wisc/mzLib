// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaModification.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;

namespace Proteomics
{
    public class ChemicalFormulaModification : Modification, IHasChemicalFormula
    {
        /// <summary>
        /// The Chemical Formula of this modifications
        /// </summary>
        public ChemicalFormula ThisChemicalFormula { get; private set; }

        public ChemicalFormulaModification(ChemicalFormula chemicalFormula)
            : this(chemicalFormula, ModificationSites.Any)
        {
        }

        public ChemicalFormulaModification(ChemicalFormula chemicalFormula, ModificationSites sites)
            : this(chemicalFormula, "", sites)
        {
            Name = ThisChemicalFormula.Formula;
        }

        public ChemicalFormulaModification(ChemicalFormula chemicalFormula, string name)
            : this(chemicalFormula, name, ModificationSites.Any)
        {
        }

        public ChemicalFormulaModification(ChemicalFormula chemicalFormula, string name, ModificationSites sites)
            : base(chemicalFormula.MonoisotopicMass, name, sites)
        {
            ThisChemicalFormula = chemicalFormula;
        }

        public ChemicalFormulaModification(ChemicalFormulaModification other)
            : this(new ChemicalFormula(other.ThisChemicalFormula), other.Name, other.Sites)
        {
        }
    }
}