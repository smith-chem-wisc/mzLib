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
using MzLibUtil;
using System;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace MzIdentML
{
    public class MzidIdentifications : IIdentifications
    {

        #region Private Fields

        private readonly mzIdentML.Generated.MzIdentMLType dd;
        private readonly mzIdentML110.Generated.MzIdentMLType dd110;

        #endregion Private Fields

        #region Public Constructors

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

        #endregion Public Constructors

        #region Public Properties

        public Tolerance ParentTolerance
        {
            get
            {
                try
                {
                    var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value)) :
                           new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
                catch
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value)) :
                           new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
            }
        }

        public Tolerance FragmentTolerance
        {
            get
            {
                try
                {
                    var hm = dd.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value)) :
                           new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
                }
                catch
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           new Tolerance(ToleranceUnit.Absolute, Convert.ToDouble(hm[0].value)) :
                           new Tolerance(ToleranceUnit.PPM, Convert.ToDouble(hm[0].value));
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

        #endregion Public Properties

        #region Public Methods

        public double CalculatedMassToCharge(int sirIndex)
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

        public int ChargeState(int sirIndex)
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

        public double ExperimentalMassToCharge(int sirIndex)
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

        public bool IsDecoy(int sirIndex)
        {
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        return ok.isDecoy;
                    }
                }
                return false;
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        return ok.isDecoy;
                    }
                }
                return false;
            }
        }

        public double QValue(int sirIndex)
        {
            try
            {
                var cvParam = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].cvParam.Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return cvParam == null? -1 : Convert.ToDouble(cvParam.value);
            }
            catch
            {
                var cvParam = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].cvParam.Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return cvParam == null ? -1 : Convert.ToDouble(cvParam.value);
            }
        }

        public string ModificationAcession(int sirIndex, int i)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.Modification[i].cvParam[0].accession;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.Modification[i].cvParam[0].accession;
                                break;
                            }
                        }
                    }
                }
            }
            return s;
        }

        public string ModificationDictionary(int sirIndex, int i)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.Modification[i].cvParam[0].cvRef;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.Modification[i].cvParam[0].cvRef;
                                break;
                            }
                        }
                    }
                }
            }
            return s;
        }

        public int ModificationLocation(int sirIndex, int i)
        {
            int modLoc = -1;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                modLoc = ok2.Modification[i].location;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                modLoc = ok2.Modification[i].location;
                                break;
                            }
                        }
                    }
                }
            }
            return modLoc;
        }

        public int NumModifications(int sirIndex)
        {
            int numMod = 0;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                if (ok2.Modification == null)
                                    break;
                                numMod = ok2.Modification.Length;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                if (ok2.Modification == null)
                                    break;
                                numMod = ok2.Modification.Length;
                                break;
                            }
                        }
                    }
                }
            }
            return numMod;
        }

        public string PeptideSequenceWithoutModifications(int sirIndex)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.PeptideSequence;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.PeptideSequence;
                                break;
                            }
                        }
                    }
                }
            }
            return s;
        }

        public string Ms2SpectrumID(int sirIndex)
        {
            string ms2id = null;
            try
            {
                if (dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
                || dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                {
                    ms2id = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                }
                else if (dd.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                {
                    ms2id = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                }
            }
            catch
            {
                if (dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
       || dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                {
                    ms2id = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                }
                else if (dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                {
                    ms2id = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                }
            }
            return ms2id;
        }

        public float[] MatchedIons(int sirIndex, int i)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].Fragmentation[i].FragmentArray[0].values;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].Fragmentation[i].FragmentArray[0].values;
            }
        }

        public int MatchedIonCounts(int sirIndex, int i)
        {
            try
            {
                return dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].Fragmentation[i].FragmentArray[0].values.Length;
            }
            catch
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].Fragmentation[i].FragmentArray[0].values.Length;
            }
        }

        public string ProteinAccession(int sirIndex)
        {
            string s = null;

            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.DBSequence)
                        {
                            if (ok2.id.Equals(ok.dBSequence_ref))
                            {
                                s = ok2.accession;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.DBSequence)
                        {
                            if (ok2.id.Equals(ok.dBSequence_ref))
                            {
                                s = ok2.accession;
                                break;
                            }
                        }
                    }
                }
            }
            return s;
        }

        public string ProteinFullName(int sirIndex)
        {
            string s = null;

            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd.SequenceCollection.DBSequence)
                        {
                            if (ok2.id.Equals(ok.dBSequence_ref))
                            {
                                s = ok2.name;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.DBSequence)
                        {
                            if (ok2.id.Equals(ok.dBSequence_ref))
                            {
                                s = ok2.name;
                                break;
                            }
                        }
                    }
                }
            }
            return s;
        }

        public int StartResidueInProtein(int sirIndex)
        {
            int startResidue = -1;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        startResidue = ok.start;
                        break;
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        startResidue = ok.start;
                        break;
                    }
                }
            }
            return startResidue;
        }

        public int EndResidueInProtein(int sirIndex)
        {
            int endResidue = -1;
            try
            {
                string peptideEvidenceRef = dd.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        endResidue = ok.end;
                        break;
                    }
                }
            }
            catch
            {
                string peptideEvidenceRef = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[0].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        endResidue = ok.end;
                        break;
                    }
                }
            }
            return endResidue;
        }

        #endregion Public Methods

    }
}