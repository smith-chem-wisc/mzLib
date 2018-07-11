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

using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    public class FakeMsDataFile : MsDataFile
    {
        public FakeMsDataFile(MsDataScan[] FakeScans) : base(FakeScans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
        {
            this.Scans = FakeScans;
        }

        public override int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            int ok = Array.BinarySearch(Scans.Select(b => b.RetentionTime).ToArray(), retentionTime);
            if (ok < 0)
                ok = ~ok;
            return ok + 1;
        }

        public override IEnumerable<MsDataScan> GetMS1Scans()
        {
            throw new NotImplementedException();
        }

        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }
    }
}