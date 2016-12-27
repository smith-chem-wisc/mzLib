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
        private mzIdentML.Generated.MzIdentMLType dd = null;
        private mzIdentML110.Generated.MzIdentMLType dd110 = null;

        public MzidIdentifications(string mzidFile)
        {
            try
            {
                using (Stream stream = new FileStream(mzidFile, FileMode.Open))
                {
                    XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML.Generated.MzIdentMLType));
                    // Read the XML file into the variable
                    dd = _indexedSerializer.Deserialize(stream) as mzIdentML.Generated.MzIdentMLType;
                }
            }
            catch
            {
                using (Stream stream = new FileStream(mzidFile, FileMode.Open))
                {
                    XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType));
                    // Read the XML file into the variable
                    dd110 = _indexedSerializer.Deserialize(stream) as mzIdentML110.Generated.MzIdentMLType;
                }
            }
        }

        public double calculatedMassToCharge(int sirIndex)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].calculatedMassToCharge;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].calculatedMassToCharge;
            }
        }

        public int chargeState(int sirIndex)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].chargeState;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].chargeState;
            }
        }

        public double experimentalMassToCharge(int sirIndex)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].experimentalMassToCharge;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].experimentalMassToCharge;
            }
        }

        public Tolerance parentTolerance
        {
            get
            {
                try
                {
                    var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    if (hm[0].unitName.Equals("dalton"))
                        return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                    else
                        return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
                catch
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    if (hm[0].unitName.Equals("dalton"))
                        return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                    else
                        return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
            }
        }

        public Tolerance fragmentTolerance
        {
            get
            {
                try
                {
                    var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    if (hm[0].unitName.Equals("dalton"))
                        return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                    else
                        return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
                catch
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    if (hm[0].unitName.Equals("dalton"))
                        return new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value));
                    else
                        return new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
            }
        }

        public int Count
        {
            get
            {
                try
                {
                    return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
                catch
                {
                    return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
            }
        }

        public bool isDecoy(int sirIndex)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                return dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].isDecoy;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                return dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].isDecoy;
            }
        }

        public bool passThreshold(int sirIndex)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].passThreshold;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].passThreshold;
            }
        }

        public string modificationAcession(int sirIndex, int i)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;
                foreach (var ok in dd.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.Modification[i].cvParam[0].accession;
                }
                return null;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;
                foreach (var ok in dd110.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.Modification[i].cvParam[0].accession;
                }
                return null;
            }
        }

        public string modificationDictionary(int sirIndex, int i)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.Modification[i].cvParam[0].cvRef;
                }

                return null;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd110.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.Modification[i].cvParam[0].cvRef;
                }

                return null;
            }
        }

        public int modificationLocation(int sirIndex, int i)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))

                        return ok.Modification[i].location;
                }

                return -1;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd110.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))

                        return ok.Modification[i].location;
                }

                return -1;
            }
        }

        public int NumModifications(int sirIndex)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                    {
                        if (ok.Modification == null)
                            return 0;
                        return ok.Modification.Length;
                    }
                }
                return -1;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd110.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                    {
                        if (ok.Modification == null)
                            return 0;
                        return ok.Modification.Length;
                    }
                }
                return -1;
            }
        }

        public string PeptideSequenceWithoutModifications(int sirIndex)
        {
            try
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.PeptideSequence;
                }

                return null;
            }
            catch
            {
                int peptideEvidenceIndex = GetLastNumberFromString(dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref);
                var peptideRef = dd110.SequenceCollection.PeptideEvidence[peptideEvidenceIndex - 1].peptide_ref;

                foreach (var ok in dd110.SequenceCollection.Peptide)
                {
                    if (ok.id.Equals(peptideRef))
                        return ok.PeptideSequence;
                }

                return null;
            }
        }

        public string ms2spectrumID(int sirIndex)
        {
            try
            {
                if (dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
                || dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                {
                    return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                }
                else if (dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                {
                    return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                }
                else
                    return null;
            }
            catch
            {
                if (dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
       || dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                {
                    return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                }
                else if (dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                {
                    return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                }
                else
                    return null;
            }
        }

        private static int GetLastNumberFromString(string s)
        {
            return Convert.ToInt32(Regex.Match(s, @"\d+$").Value);
        }
    }
}