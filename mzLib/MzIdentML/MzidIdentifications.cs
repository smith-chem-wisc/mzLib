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
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace MzIdentML
{
    public class MzidIdentifications : IIdentifications
    {

        private readonly mzIdentML110.Generated.MzIdentMLType110 dd110;
        private readonly mzIdentML111.Generated.MzIdentMLType111 dd111;
        private readonly mzIdentML120.Generated.MzIdentMLType120 dd120;



        public MzidIdentifications(string mzidFile)
        {
            try
            {
                using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
                    // Read the XML file into the variable
                    dd110 = _indexedSerializer.Deserialize(stream) as mzIdentML110.Generated.MzIdentMLType110;
                }
            }
            catch
            {
                try
                {
                    using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
                        // Read the XML file into the variable
                        dd111 = _indexedSerializer.Deserialize(stream) as mzIdentML111.Generated.MzIdentMLType111;
                    }
                }
                catch
                {
                    using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
                        // Read the XML file into the variable
                        dd120 = _indexedSerializer.Deserialize(stream) as mzIdentML120.Generated.MzIdentMLType120;
                    }
                }

                
            }
        }



        public Tolerance ParentTolerance
        {
            get
            {
                try
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                catch
                {
                    try
                    {
                        var hm = dd111.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                        return hm[0].unitName.Equals("dalton") ?
                               (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                               new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                    }
                    catch
                    {
                        var hm = dd120.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                        return hm[0].unitName.Equals("dalton") ?
                               (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                               new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                    }
                }
            }
        }

        public Tolerance FragmentTolerance
        {
            get
            {
                try
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                catch
                {
                    try
                    {
                        var hm = dd111.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                        return hm[0].unitName.Equals("dalton") ?
                               (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                               new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                    }
                    catch
                    {
                        var hm = dd120.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                        return hm[0].unitName.Equals("dalton") ?
                               (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                               new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                    }
                }
            }
        }

        public int Count
        {
            get
            {
                try
                {
                    return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
                catch
                {
                    try
                    {
                        return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                    }
                    catch
                    {
                        return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                    }
                }
            }
        }


        #region Public Methods

        public double CalculatedMassToCharge(int sirIndex, int siiIndex)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
                }
            }
        }

        public int ChargeState(int sirIndex, int siiIndex)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
                }
            }
        }

        public double ExperimentalMassToCharge(int sirIndex, int siiIndex)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
                }
            }
        }

        public bool IsDecoy(int sirIndex, int siiIndex)
        {
            //if any of the peptide evidences is decoy, is decoy
            try
            {
                foreach (mzIdentML110.Generated.PeptideEvidenceRefType pe in dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            if (!ok.isDecoy) return false;
                        }
                    }
                }
                return true;
            }
            catch
            {
                try
                {

                    foreach (mzIdentML111.Generated.PeptideEvidenceRefType pe 
                        in dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                        {
                            //if (ok.id.Equals(peptideEvidenceRef))
                            //{
                            //    if (!ok.isDecoy) return false;
                            //}

                            if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                            {
                                return false;
                            }
                        }
                    }
                    return true;
                }
                catch
                {
                    foreach (mzIdentML120.Generated.PeptideEvidenceRefType pe 
                        in dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                        {
                            //if (ok.id.Equals(peptideEvidenceRef))
                            //{
                            //    if (!ok.isDecoy) return false;
                            //}

                            if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                            {
                                return false;
                            }

                        }
                    }
                    return true;
                }
            }
        }

        public double QValue(int sirIndex, int siiIndex)
        {
            try
            {
                var cvParam = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam.
                    Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return cvParam == null ? -1 : Convert.ToDouble(cvParam.value, CultureInfo.InvariantCulture);
            }
            catch
            {
                try
                {
                    var cvParam = dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam.
                        Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                    return cvParam == null ? -1 : Convert.ToDouble(cvParam.value, CultureInfo.InvariantCulture);
                }
                catch
                {
                    var cvParam = dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam.
                        Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                    return cvParam == null ? -1 : Convert.ToDouble(cvParam.value, CultureInfo.InvariantCulture);
                }

            }
        }

        public int NumPSMsFromScan(int sirIndex)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
                }
            }
        }

        public string ModificationAcession(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
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

                
            }
            return s;
        }

        public string ModificationValue(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                s = ok2.Modification[i].cvParam[0].value;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
                            {
                                if (ok2.id.Equals(ok.peptide_ref))
                                {
                                    s = ok2.Modification[i].cvParam[0].value;
                                    break;
                                }
                            }
                        }
                    }
                }
                catch
                {
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
                            {
                                if (ok2.id.Equals(ok.peptide_ref))
                                {
                                    s = ok2.Modification[i].cvParam[0].value;
                                    break;
                                }
                            }
                        }
                    }
                }
                
            }
            return s;
        }

        public string ModificationDictionary(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
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
                
            }
            return s;
        }

        public int ModificationLocation(int sirIndex, int siiIndex, int i)
        {
            int modLoc = -1;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
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
                
            }
            return modLoc;
        }

        public double ModificationMass(int sirIndex, int siiIndex, int i)
        {
            double modMass = -1;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd110.SequenceCollection.Peptide)
                        {
                            if (ok2.id.Equals(ok.peptide_ref))
                            {
                                modMass = ok2.Modification[i].monoisotopicMassDelta;
                                break;
                            }
                        }
                    }
                }
            }
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
                            {
                                if (ok2.id.Equals(ok.peptide_ref))
                                {
                                    modMass = ok2.Modification[i].monoisotopicMassDelta;
                                    break;
                                }
                            }
                        }
                    }
                }
                catch
                {
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
                            {
                                if (ok2.id.Equals(ok.peptide_ref))
                                {
                                    modMass = ok2.Modification[i].monoisotopicMassDelta;
                                    break;
                                }
                            }
                        }
                    }
                }
                
            }
            return modMass;
        }

        public int NumModifications(int sirIndex, int siiIndex)
        {
            int numMod = 0;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
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
                
            }
            return numMod;
        }

        public string PeptideSequenceWithoutModifications(int sirIndex, int siiIndex)
        {
            string s = null;
            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.Peptide)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.Peptide)
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
                
            }
            return s;
        }

        public string Ms2SpectrumID(int sirIndex)
        {
            string ms2id = null;
            try
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
            catch
            {
                try
                {
                    if (dd111.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
       || dd111.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                    {
                        ms2id = dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                    }
                    else if (dd111.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                    {
                        ms2id = dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                    }
                }
                catch
                {
                    if (dd120.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Thermo RAW format")
       || dd120.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("mzML format"))
                    {
                        ms2id = dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].spectrumID;
                    }
                    else if (dd120.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam.name.Equals("Mascot MGF format"))
                    {
                        ms2id = dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].cvParam[0].value;
                    }
                }

                
            }
            return ms2id;
        }

        public float[] MatchedIons(int sirIndex, int siiIndex, int i)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
                }
            }
        }

        public int MatchedIonCounts(int sirIndex, int siiIndex, int i)
        {
            try
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
            }
            catch
            {
                try
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
                }
                catch
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
                }
            }
        }

        public string ProteinAccession(int sirIndex, int siiIndex)
        {
            string s = null;

            try
            {
                string peptideEvidenceRef = 
                    dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
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
            catch
            {
                try
                {
                    string peptideEvidenceRef = 
                        dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd111.SequenceCollection.DBSequence)
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
                    string peptideEvidenceRef = 
                        dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd120.SequenceCollection.DBSequence)
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
                
            }
            return s;
        }

        public string ProteinFullName(int sirIndex, int siiIndex)
        {
            string s = "";

            try
            {
                foreach (mzIdentML110.Generated.PeptideEvidenceRefType pe 
                    in dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd110.SequenceCollection.DBSequence)
                            {
                                if (ok2.id.Equals(ok.dBSequence_ref))
                                {
                                    if (s.Length != 0) s += " or ";
                                    s += ok2.name;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            catch
            {
                try
                {
                    foreach (mzIdentML111.Generated.PeptideEvidenceRefType pe 
                        in dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                foreach (var ok2 in dd111.SequenceCollection.DBSequence)
                                {
                                    if (ok2.id.Equals(ok.dBSequence_ref))
                                    {
                                        if (s.Length != 0) s += " or ";
                                        s += ok2.name;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                catch
                {
                    foreach (mzIdentML120.Generated.PeptideEvidenceRefType pe 
                        in dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                foreach (var ok2 in dd120.SequenceCollection.DBSequence)
                                {
                                    if (ok2.id.Equals(ok.dBSequence_ref))
                                    {
                                        if (s.Length != 0) s += " or ";
                                        s += ok2.name;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
            }
            return s;
        }

        public string StartResidueInProtein(int sirIndex, int siiIndex)
        {
            string startResidue = "";
            try
            {
                foreach (mzIdentML110.Generated.PeptideEvidenceRefType pe 
                    in dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            if (startResidue.Length != 0) startResidue += " or ";
                            startResidue += ok.start;
                            break;
                        }
                    }
                }
            }
            catch
            {
                try
                {

                    foreach (mzIdentML111.Generated.PeptideEvidenceRefType pe 
                        in dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                if (startResidue.Length != 0) startResidue += " or ";
                                startResidue += ok.start;
                                break;
                            }
                        }
                    }
                }
                catch
                {
                    foreach (mzIdentML120.Generated.PeptideEvidenceRefType pe 
                        in dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                if (startResidue.Length != 0) startResidue += " or ";
                                startResidue += ok.start;
                                break;
                            }
                        }
                    }
                }
                
            }
            return startResidue;
        }

        public string EndResidueInProtein(int sirIndex, int siiIndex)
        {
            string endResidue = "";
            try
            {
                foreach (mzIdentML110.Generated.PeptideEvidenceRefType pe 
                    in dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            if (endResidue.Length != 0) endResidue += " or ";
                            endResidue += ok.end;
                            break;
                        }
                    }
                }
            }
            catch
            {
                try
                {
                    foreach (mzIdentML111.Generated.PeptideEvidenceRefType pe 
                        in dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                if (endResidue.Length != 0) endResidue += " or ";
                                endResidue += ok.end;
                                break;
                            }
                        }
                    }
                }
                catch
                {
                    foreach (mzIdentML120.Generated.PeptideEvidenceRefType pe 
                        in dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                    {
                        string peptideEvidenceRef = pe.peptideEvidence_ref;
                        foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                        {
                            if (ok.id.Equals(peptideEvidenceRef))
                            {
                                if (endResidue.Length != 0) endResidue += " or ";
                                endResidue += ok.end;
                                break;
                            }
                        }
                    }
                }
                
            }
            return endResidue;
        }

        #endregion Public Methods
    }
}