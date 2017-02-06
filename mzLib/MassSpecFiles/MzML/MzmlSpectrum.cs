// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (DefaultMzSpectrum.cs) is part of MassSpectrometry.
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

using System;
using MassSpectrometry;
using Spectra;

namespace IO.MzML
{
    public class MzmlMzSpectrum : MzSpectrum<MzmlPeak>
    {
        public MzmlMzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        public override Spectrum<MzmlPeak> CreateSpectrumFromTwoArrays(double[] item1, double[] item2, bool v)
        {
            return GetMzSpectrumFromTwoArrays(item1, item2, v);
        }

        public override MzSpectrum<MzmlPeak> GetMzSpectrumFromTwoArrays(double[] item1, double[] item2, bool v)
        {
            return new MzmlMzSpectrum(item1, item2, v);
        }


        public new MzmlMzSpectrum NewSpectrumApplyFunctionToX(Func<double, double> convertor)
        {
            var ok = ApplyFunctionToX(convertor);
            return new MzmlMzSpectrum(ok.Item1, ok.Item2, false);
        }
    }
}