// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (IMzSpectrum.cs) is part of MassSpectrometry.
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

using MzLibUtil;
using Spectra;
using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public interface IMzSpectrum<out TPeak> : ISpectrum<TPeak>
        where TPeak : IMzPeak
    {
        #region Public Properties

        new MzRange Range { get; }

        #endregion Public Properties

        #region Public Methods

        void ReplaceXbyApplyingFunction(Func<IMzPeak, double> convertor);

        byte[] Get64BitXarray();

        byte[] Get64BitYarray();

        IEnumerable<Tuple<List<IMzPeak>, int>> Deconvolute(MzRange range, int maxAssumedChargeState, Tolerance massTolerance, double intensityRatio);

        #endregion Public Methods
    }
}