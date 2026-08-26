// Copyright 2017 Stefan Solntsev
//
// This file (SourceFile.cs) is part of MassSpectrometry.
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
using System.IO;
using MzLibUtil;

namespace MassSpectrometry
{
    public class SourceFile
    {
        /// <summary>
        /// Written to mzML's <c>sourceFile/@location</c> when no usable <see cref="Uri"/> is available,
        /// and used by the mzML reader for the same case. The attribute is <c>use="required"</c> in the
        /// schema, so it cannot simply be omitted. Shared so the two sides cannot drift apart.
        /// </summary>
        public const string UnknownLocation = "file:///unknown-source";

        /// <summary>
        /// Written to mzML's <c>sourceFile/@name</c> when <see cref="FileName"/> is null. Also
        /// <c>use="required"</c>, and omitting it produced a schema-invalid file.
        /// </summary>
        public const string UnknownName = "unknown-source";

        public SourceFile(string nativeIdFormat, string massSpectrometerFileFormat, string checkSum, string fileChecksumType, string id)
        {
            NativeIdFormat = nativeIdFormat;
            MassSpectrometerFileFormat = massSpectrometerFileFormat;
            CheckSum = checkSum;
            FileChecksumType = fileChecksumType;
            Id = id;
        }

        public SourceFile(string nativeIdFormat, string massSpectrometerFileFormat, string checkSum, string fileChecksumType, string filePath, string id)
        : this(nativeIdFormat, massSpectrometerFileFormat, checkSum, fileChecksumType, id)
        {
            Uri.TryCreate(Directory.GetParent(filePath).FullName, UriKind.Absolute, out Uri result);
            this.Uri = result;
            this.FileName = Path.GetFileName(filePath);
        }

        public SourceFile(string nativeIdFormat, string massSpectrometerFileFormat, string checkSum, string fileChecksumType, Uri uri, string id, string fileName)
            : this(nativeIdFormat, massSpectrometerFileFormat, checkSum, fileChecksumType, id)
        {
            this.Uri = uri;
            this.FileName = fileName;
        }

        public string NativeIdFormat { get; }
        public string MassSpectrometerFileFormat { get; }
        public string CheckSum { get; }
        public string FileChecksumType { get; }
        public Uri Uri { get; }

        public string FileName { get; }
        public string Id { get; }

        /// <summary>
        /// The mass spectrometer that produced this file, e.g. NT=Q Exactive;AC=MS:1001911.
        /// Null when the source format does not record it or the reader could not determine it.
        ///
        /// Init-only, so every existing construction site is unaffected: this is additive.
        ///
        /// The two readers that can populate it do so with different completeness, and callers need
        /// to know which they have:
        ///
        ///   mzML   carries the instrument as a PSI-MS cvParam in instrumentConfigurationList, so
        ///          both <see cref="MzLibUtil.CvParam.Accession"/> and Name are set. Nothing has to
        ///          be looked up or guessed.
        ///   Thermo RAW records only a model NAME ("Orbitrap Fusion Lumos"), so Name is set and
        ///          Accession is empty. Resolving that name to an accession needs a controlled-
        ///          vocabulary lookup, which is a separate concern from reading the file.
        ///
        /// Consumers must therefore check Accession for emptiness rather than assuming a populated
        /// InstrumentModel is fully accessioned. Per the rule on <see cref="MzLibUtil.CvParam"/>,
        /// match on Accession and never on Name.
        /// </summary>
        public CvParam InstrumentModel { get; init; }
    }
}