// Copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaTerminus.cs) is part of Proteomics.
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
    public class ChemicalFormulaTerminus : IHasChemicalFormula
    {
        public ChemicalFormulaTerminus(ChemicalFormula chemicalFormula)
        {
            ThisChemicalFormula = chemicalFormula;
        }

        public double MonoisotopicMass
        {
            get
            {
                return ThisChemicalFormula.MonoisotopicMass;
            }
        }

        public ChemicalFormula ThisChemicalFormula
        {
            get; private set;
        }
    }
}