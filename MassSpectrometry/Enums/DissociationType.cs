// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (DissociationType.cs) is part of MassSpectrometry.
//
// MassSpectrometry is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry. If not, see <http://www.gnu.org/licenses/>.

namespace MassSpectrometry
{
    public enum DissociationType
    {
        Unknown = -1,

        // The values below are identical to thermo names
        CID = 0, // MS:1000133 collision-induced dissociation

        IRMPD = 1, // MS:1000435 photodissociation
        ECD = 2, // MS:1000250 electron capture dissociation
        PQD = 3, // MS:1000599 pulsed q dissociation
        ETD = 4, // MS:1000598 electron transfer dissociation
        HCD = 5, // MS:1002481 higher energy beam-type collision-induced dissociation

        AnyActivationType = 6,

        EThcD = 7, // MS:1002631 Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)
        Custom = 8,

        ISCID = 9,
        //NPTR = 10,
        // The values above are identical to thermo names

        //ISCID = 11
    }
}