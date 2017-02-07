// Copyright 2016 Stefan Solntsev
//
// This file (FakeMsDataFile.cs) is part of MassSpectrometry.
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

using IO.MzML;
using MassSpectrometry;
using System;
using System.Linq;

namespace Test
{
    public class FakeMsDataFile : MsDataFile<IMzmlScan, MzmlMzSpectrum, MzmlPeak>
    {
        #region Public Constructors

        public FakeMsDataFile(string filePath, IMzmlScan[] FakeScans) : base(filePath, MsDataFileType.UnKnown)
        {
            this.Scans = FakeScans;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            int ok = Array.BinarySearch(Scans.Select(b => b.RetentionTime).ToArray(), retentionTime);
            if (ok < 0)
                ok = ~ok;
            return ok + 1;
        }

        public override void Open()
        {
        }

        public override void Close()
        {
        }

        #endregion Public Methods

        #region Protected Methods

        protected override int GetNumSpectra()
        {
            return Scans.Count();
        }

        protected override IMzmlScan GetMsDataOneBasedScanFromFile(int oneBasedSpectrumNumber)
        {
            return Scans[oneBasedSpectrumNumber - 1];
        }

        #endregion Protected Methods
    }
}