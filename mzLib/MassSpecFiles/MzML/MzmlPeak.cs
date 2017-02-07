// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MZPeak.cs) is part of MassSpectrometry.
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

using MassSpectrometry;

namespace IO.MzML
{
    /// <summary>
    /// A peak in a mass spectrum that has a well defined m/z and intenisty value
    /// </summary>
    public class MzmlPeak : MzPeak
    {

        #region Public Constructors

        public MzmlPeak(double mz, double intensity) : base(mz, intensity)
        {
        }

        #endregion Public Constructors

    }
}