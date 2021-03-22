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
    /// <summary>
    /// See https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo for PSI-MS list of dissociation types.
    /// search for "is_a: MS:1000044 ! dissociation method" and "is_a: MS:1000422 ! beam-type collision-induced dissociation"
    /// </summary>
    public enum DissociationType
    {
        /// <summary>
        /// id: MS:1000133 collision-induced dissociation
        /// </summary>
        CID,

        /// <summary>
        /// id: MS:1000134 plasma desorption
        /// </summary>
        PD,

        /// <summary>
        /// id: MS:1000135 post-source decay
        /// </summary>
        PSD,

        /// <summary>
        /// id: MS:1000136 surface-induced dissociation
        /// </summary>
        SID,

        /// <summary>
        /// id: MS:1000242 blackbody infrared radiative dissociation
        /// </summary>
        BIRD,

        /// <summary>
        /// MS:1000250 electron capture dissociation
        /// </summary>
        ECD,

        /// <summary>
        /// MS:1000262 infrared multiphoton dissociation
        /// </summary>
        IRMPD,

        /// <summary>
        /// id: MS:1000282 sustained off-resonance irradiation
        /// </summary>
        SORI,

        /// <summary>
        /// MS:1000435 photodissociation
        /// </summary>
        MPD,

        /// <summary>
        /// MS:1000598 electron transfer dissociation
        /// </summary>
        ETD,

        /// <summary>
        /// MS:1000599 pulsed q dissociation
        /// </summary>
        PQD,

        /// <summary>
        /// id: MS:1001880 in-source collision-induced dissociation
        /// </summary>
        ISCID,

        /// <summary>
        /// MS:1000422 beam-type collision-induced dissociation
        /// </summary>
        HCD,

        /// <summary>
        /// MS:1002631 Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)
        /// </summary>
        EThcD,

        /// <summary>
        /// ultraviolet photodissociation
        /// </summary>
        UVPD,

        /// <summary>
        /// negative-mode electron transfer dissociation
        /// </summary>
        NETD,

        /// <summary>
        /// low-resolution (ion-trap) CID
        /// </summary>
        LowCID,

        Unknown,
        AnyActivationType,
        Custom,
        
        /// <summary>
        /// Placeholder used to communicate to MetaMorpheus that the dissociation type 
        /// should be taken from the scan header instead of using a fixed, constant dissociation type
        /// </summary>
        Autodetect
    }
}