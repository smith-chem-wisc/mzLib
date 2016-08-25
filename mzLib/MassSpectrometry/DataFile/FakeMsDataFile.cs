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

using Spectra;
using System;
using System.Linq;

namespace MassSpectrometry
{
    public class FakeMsDataFile : MsDataFile<IMzSpectrum<MzPeak>>
    {
        private MsDataScan<IMzSpectrum<MzPeak>>[] FakeScans;
        public FakeMsDataFile(string filePath, MsDataScan<IMzSpectrum<MzPeak>>[] FakeScans) : base(filePath, true, MsDataFileType.UnKnown)
        {
            this.FakeScans = FakeScans;
        }

        public override int GetSpectrumNumber(double retentionTime)
        {
            int ok = Array.BinarySearch(Scans.Select(b => b.RetentionTime).ToArray(), retentionTime);
            if (ok < 0)
                ok = ~ok;
            return ok + FirstSpectrumNumber;
        }

        public override void Open()
        {
        }

        protected override int GetFirstSpectrumNumber()
        {
            return 1;
        }

        protected override int GetLastSpectrumNumber()
        {
            return FakeScans.Count();
        }

        protected override MsDataScan<IMzSpectrum<MzPeak>> GetMsDataScanFromFile(int spectrumNumber)
        {
            return FakeScans[spectrumNumber - 1];
        }
    }
}
