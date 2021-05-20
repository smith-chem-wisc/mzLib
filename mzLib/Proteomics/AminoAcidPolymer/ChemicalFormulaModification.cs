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

namespace Proteomics.AminoAcidPolymer
{
    public class OldSchoolChemicalFormulaModification : OldSchoolModification, IHasChemicalFormula
    {
        public OldSchoolChemicalFormulaModification(ChemicalFormula chemicalFormula)
            : this(chemicalFormula, ModificationSites.Any)
        {
        }

        public OldSchoolChemicalFormulaModification(ChemicalFormula chemicalFormula, ModificationSites sites)
            : this(chemicalFormula, "", sites)
        {
            Name = ThisChemicalFormula.Formula;
        }

        public OldSchoolChemicalFormulaModification(ChemicalFormula chemicalFormula, string name)
            : this(chemicalFormula, name, ModificationSites.Any)
        {
        }

        public OldSchoolChemicalFormulaModification(ChemicalFormula chemicalFormula, string name, ModificationSites sites)
            : base(chemicalFormula.MonoisotopicMass, name, sites)
        {
            ThisChemicalFormula = chemicalFormula;
        }

        public OldSchoolChemicalFormulaModification(OldSchoolChemicalFormulaModification other)
            : this(ChemicalFormula.ParseFormula(other.ThisChemicalFormula.Formula), other.Name, other.Sites)
        {
        }

        /// <summary>
        /// The Chemical Formula of this modifications
        /// </summary>
        public ChemicalFormula ThisChemicalFormula { get; private set; }
    }
}