// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (IMsDataFile.cs) is part of MassSpectrometry.
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

using System.IO;

namespace MassSpectrometry
{
    public class SourceFile
    {
        #region Public Constructors

        public SourceFile()
        {
            this.NativeIdFormat = @"no nativeID format";
            this.MassSpectrometerFileFormat = @"mzML format";
            this.CheckSum = @"";
            this.FileChecksumType = @"SHA-1";
            this.FilePath = @"C:\undefined.mzML";
            this.Id = @"undefined.mzML";
        }

        public SourceFile(string nativeIdFormat, string massSpectrometerFileFormat, string checkSum, string fileChecksumType, string filePath, string id)
        {
            this.NativeIdFormat = nativeIdFormat;
            this.MassSpectrometerFileFormat = massSpectrometerFileFormat;
            this.CheckSum = checkSum;
            this.FileChecksumType = fileChecksumType;
            this.FilePath = filePath;
            this.Id = id;
        }

        #endregion Public Constructors

        #region Public Properties

        public string NativeIdFormat { get; }
        public string MassSpectrometerFileFormat { get; }
        public string CheckSum { get; }
        public string FileChecksumType { get; }
        public string FileLocation { get { return Directory.GetParent(FilePath).FullName; } }
        public string FileName { get { return Path.GetFileName(FilePath); } }
        public string FilePath { get; }
        public string Id { get; }

        #endregion Public Properties
    }
}