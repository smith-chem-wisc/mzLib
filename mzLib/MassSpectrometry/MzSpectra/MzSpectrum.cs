// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MzSpectrum.cs) is part of MassSpectrometry.
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
using System.IO;

namespace MassSpectrometry
{
    public abstract class MzSpectrum<TPeak> : Spectrum<TPeak>, IMzSpectrum<TPeak>
        where TPeak : IMzPeak
    {

        #region Protected Constructors

        protected MzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        protected MzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Protected Constructors

        #region Public Properties

        new public MzRange Range
        {
            get
            {
                return new MzRange(FirstX, LastX);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static byte[] Get64Bitarray(IEnumerable<double> yArray)
        {
            var mem = new MemoryStream();
            foreach (var okk in yArray)
            {
                byte[] ok = BitConverter.GetBytes(okk);
                mem.Write(ok, 0, ok.Length);
            }
            mem.Position = 0;
            return mem.ToArray();
        }

        public byte[] Get64BitYarray()
        {
            return Get64Bitarray(YArray);
        }

        public byte[] Get64BitXarray()
        {
            return Get64Bitarray(XArray);
        }

        public void ReplaceXbyApplyingFunction(Func<IMzPeak, double> convertor)
        {
            for (int i = 0; i < Size; i++)
                XArray[i] = convertor(this[i]);
            peakWithHighestY = default(TPeak);
            peakList = new TPeak[Size];
        }

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Size);
        }

        #endregion Public Methods

    }
}