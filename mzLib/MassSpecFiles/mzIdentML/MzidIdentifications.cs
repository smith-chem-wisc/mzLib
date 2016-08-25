// Copyright 2016 Stefan Solntsev
// 
// This file (MzidIdentifications.cs) is part of MassSpecFiles.
// 
// MassSpecFiles is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// MassSpecFiles is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpecFiles. If not, see <http://www.gnu.org/licenses/>.

using MassSpectrometry;
using Spectra;
using System;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml.Serialization;

namespace MzIdentML
{
    public class MzidIdentifications : Identifications
    {
        private mzIdentML.Generated.MzIdentMLType dd;
        public MzidIdentifications(string mzidFile)
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML.Generated.MzIdentMLType));
            Stream stream = new FileStream(mzidFile, FileMode.Open);
            // Read the XML file into the variable
            dd = _indexedSerializer.Deserialize(stream) as mzIdentML.Generated.MzIdentMLType;
        }

        public double calculatedMassToCharge(int matchIndex)
        {
            return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[matchIndex].SpectrumIdentificationItem[0].calculatedMassToCharge;
        }

        public int chargeState(int matchIndex)
        {
            return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[matchIndex].SpectrumIdentificationItem[0].chargeState;
        }

        public double experimentalMassToCharge(int matchIndex)
        {
            return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[matchIndex].SpectrumIdentificationItem[0].experimentalMassToCharge;
        }

        public Tolerance parentTolerance
        {
            get
            {
                var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                if (hm[0].unitName.Equals("dalton"))
                    return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                else
                    return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
            }
        }
        public Tolerance fragmentTolerance
        {
            get
            {
                var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                if (hm[0].unitName.Equals("dalton"))
                    return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                else
                    return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
            }
        }
        public int Count
        {
            get
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
            }
        }

        public bool isDecoy(int matchIndex)
        {
            return dd.SequenceCollection.PeptideEvidence[matchIndex].isDecoy;
        }

        public string modificationAcession(int matchIndex, int i)
        {
            return dd.SequenceCollection.Peptide[matchIndex].Modification[i].cvParam[0].accession;
        }

        public string modificationDictionary(int matchIndex, int i)
        {
            return dd.SequenceCollection.Peptide[matchIndex].Modification[i].cvParam[0].cvRef;
        }

        public int modificationLocation(int matchIndex, int i)
        {
            return dd.SequenceCollection.Peptide[matchIndex].Modification[i].location;
        }

        public int NumModifications(int matchIndex)
        {
            if (dd.SequenceCollection.Peptide[matchIndex].Modification == null)
                return 0;
            return dd.SequenceCollection.Peptide[matchIndex].Modification.Length;
        }

        public string PeptideSequenceWithoutModifications(int matchIndex)
        {
            return dd.SequenceCollection.Peptide[matchIndex].PeptideSequence;
        }

        public int ms2spectrumIndex(int matchIndex)
        {
            string ms2spectrumID = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[matchIndex].spectrumID;
            return GetLastNumberFromString(ms2spectrumID);
        }

        private static int GetLastNumberFromString(string s)
        {
            return Convert.ToInt32(Regex.Match(s, @"\d+$").Value);
        }

    }
}