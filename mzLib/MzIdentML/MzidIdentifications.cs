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
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml;
using System.Xml.Serialization;

namespace MzIdentML
{
    public class MzidIdentifications : IIdentifications
    {
        /// <summary>
        /// Reports the legacy ".../mzIdentML/1.1.0" namespace as the schema's ".../mzIdentML/1.1", so a
        /// document mzLib wrote before that was corrected deserializes with the generated 1.1.0 types.
        /// </summary>
        /// <remarks>
        /// An <see cref="XmlRootAttribute"/> override is NOT sufficient and is actively worse than
        /// failing: it renames only the root, so the root deserializes while every child element stays
        /// bound to ".../1.1" and is silently dropped. Measured on a legacy document -- id came through,
        /// cvList and AnalysisSoftwareList came back null. Remapping in the reader fixes all depths.
        /// A compliant ".../1.1" document passes through untouched, so this is safe either way.
        /// </remarks>
        private sealed class LegacyMzidNamespaceReader : XmlTextReader
        {
            private const string LegacyNamespace = "http://psidev.info/psi/pi/mzIdentML/1.1.0";
            private const string SchemaNamespace = "http://psidev.info/psi/pi/mzIdentML/1.1";

            // XmlTextReader defaults to DtdProcessing.Parse, whereas the reader XmlSerializer builds for
            // the Stream overload used by every other arm prohibits DTDs. Without this, a document the
            // other arms refuse for its DTD falls through to this one and is accepted -- in either
            // namespace, since the remap is a no-op for ".../1.1".
            public LegacyMzidNamespaceReader(Stream stream) : base(stream)
            {
                DtdProcessing = DtdProcessing.Prohibit;
            }

            public override string NamespaceURI =>
                base.NamespaceURI == LegacyNamespace ? SchemaNamespace : base.NamespaceURI;
        }


        private readonly mzIdentML110.Generated.MzIdentMLType110 dd110;
        private readonly mzIdentML111.Generated.MzIdentMLType111 dd111;
        private readonly mzIdentML120.Generated.MzIdentMLType120 dd120;
        private readonly mzIdentML130.Generated.MzIdentMLType130 dd130;



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
                    try
                    {
                        using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                        {
                            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
                            // Read the XML file into the variable
                            dd120 = _indexedSerializer.Deserialize(stream) as mzIdentML120.Generated.MzIdentMLType120;
                        }
                    }
                    catch
                    {
                        try
                        {
                            using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                            {
                                XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML130.Generated.MzIdentMLType130));
                                // Read the XML file into the variable
                                dd130 = _indexedSerializer.Deserialize(stream) as mzIdentML130.Generated.MzIdentMLType130;
                            }
                        }
                        catch
                        {
                            // Last arm: an mzIdentML 1.1.0 document declaring the namespace mzLib itself
                            // used to write, ".../mzIdentML/1.1.0", instead of the schema's
                            // ".../mzIdentML/1.1". Every .mzID this library produced before that was
                            // corrected is in the old namespace, and XmlSerializer matches namespaces
                            // exactly, so without this arm every attempt above fails and the constructor
                            // throws on files we wrote. It is last because a legacy-namespace rewrite
                            // should never pre-empt a document that parses as a real format version.
                            using (Stream stream = new FileStream(mzidFile, FileMode.Open, FileAccess.Read, FileShare.Read))
                            using (LegacyMzidNamespaceReader reader = new LegacyMzidNamespaceReader(stream))
                            {
                                XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
                                dd110 = _indexedSerializer.Deserialize(reader) as mzIdentML110.Generated.MzIdentMLType110;
                            }
                        }
                    }
                }

                
            }
        }



        public Tolerance ParentTolerance
        {
            get
            {
                if (dd110 != null)
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else if (dd111 != null)
                {
                    var hm = dd111.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else if (dd120 != null)
                {
                    var hm = dd120.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else
                {
                    var hm = dd130.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
            }
        }

        public Tolerance FragmentTolerance
        {
            get
            {
                if (dd110 != null)
                {
                    var hm = dd110.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else if (dd111 != null)
                {
                    var hm = dd111.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else if (dd120 != null)
                {
                    var hm = dd120.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
                else
                {
                    var hm = dd130.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance;
                    return hm[0].unitName.Equals("dalton") ?
                           (Tolerance)new AbsoluteTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture)) :
                           new PpmTolerance(Convert.ToDouble(hm[0].value, CultureInfo.InvariantCulture));
                }
            }
        }

        public int Count
        {
            get
            {
                if (dd110 != null)
                {
                    return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
                else if (dd111 != null)
                {
                    return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
                else if (dd120 != null)
                {
                    return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
                else
                {
                    return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult.Count();
                }
            }
        }


        #region Public Methods

        public double CalculatedMassToCharge(int sirIndex, int siiIndex)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].calculatedMassToCharge;
            }
        }

        public int ChargeState(int sirIndex, int siiIndex)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].chargeState;
            }
        }

        public double ExperimentalMassToCharge(int sirIndex, int siiIndex)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].experimentalMassToCharge;
            }
        }

        public bool IsDecoy(int sirIndex, int siiIndex)
        {
            // a match counts as a target as soon as any of its peptide evidences is a target.
            // The 1.1.0 arm expressed this as a nested if; all four now use the same predicate.
            if (dd110 != null)
            {
                foreach (mzIdentML110.Generated.PeptideEvidenceRefType pe
                    in dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd110.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
            else if (dd111 != null)
            {
                foreach (mzIdentML111.Generated.PeptideEvidenceRefType pe
                    in dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd111.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
            else if (dd120 != null)
            {
                foreach (mzIdentML120.Generated.PeptideEvidenceRefType pe
                    in dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd120.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
            else
            {
                foreach (mzIdentML130.Generated.PeptideEvidenceRefType pe
                    in dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef) && !ok.isDecoy)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
        }

        // -1 is "absent", and a q-value term with no value (the attribute is optional) is absent too.
        // Convert.ToDouble reads a null string as 0, which is indistinguishable from a perfect q-value.
        private static double QValueOf(string value) =>
            string.IsNullOrWhiteSpace(value) ? -1 : Convert.ToDouble(value, CultureInfo.InvariantCulture);

        public double QValue(int sirIndex, int siiIndex)
        {
            if (dd110 != null)
            {
                var cvParam = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam?.
                    Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return QValueOf(cvParam?.value);
            }
            else if (dd111 != null)
            {
                var cvParam = dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam?.
                    Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return QValueOf(cvParam?.value);
            }
            else if (dd120 != null)
            {
                var cvParam = dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam?.
                    Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return QValueOf(cvParam?.value);
            }
            else
            {
                var cvParam = dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].cvParam?.
                    Where(cv => cv.accession == "MS:1002354").FirstOrDefault();
                return QValueOf(cvParam?.value);
            }
        }

        public int NumPSMsFromScan(int sirIndex)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem.Count(i => i != null);
            }
        }

        public string ModificationAcession(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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

        public string ModificationValue(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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
            return s;
        }

        public string ModificationDictionary(int sirIndex, int siiIndex, int i)
        {
            string s = null;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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

        public int ModificationLocation(int sirIndex, int siiIndex, int i)
        {
            int modLoc = -1;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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

        public double ModificationMass(int sirIndex, int siiIndex, int i)
        {
            double modMass = -1;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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
            return modMass;
        }

        public int NumModifications(int sirIndex, int siiIndex)
        {
            int numMod = 0;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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

        public string PeptideSequenceWithoutModifications(int sirIndex, int siiIndex)
        {
            string s = null;
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.Peptide)
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

        // FileFormat terms are recognised by accession as well as by name. Writers carry the PSI-MS
        // accession but not always its display name -- MS-GF+ writes "mzML file" and Scaffold writes
        // "Mascot MGF file" -- and matching the name alone returned null for every result in their files.
        private const string ThermoRawFormatAccession = "MS:1000563";
        private const string MzmlFormatAccession = "MS:1000584";
        private const string MascotMgfFormatAccession = "MS:1001062";
        private const string SpectrumTitleAccession = "MS:1000796";

        // The obsolete "spectrum title", replaced_by MS:1000796, which files written against an older CV carry
        private const string ObsoleteSpectrumTitleAccession = "MS:1001416";

        private static bool IsFileFormat(string accession, string name, string expectedAccession, string expectedName) =>
            accession == expectedAccession || name == expectedName;

        /// <summary>
        /// The spectrum title among a SpectrumIdentificationResult's cvParams, found by accession. Falls back
        /// to the first cvParam, which is what was read before, only when no title term is present: Mascot
        /// Parser writes "Mascot:identity threshold" first, so reading position 0 returned the threshold.
        /// A title term that is present but has no value (the attribute is optional) is a missing title and
        /// returns null. It does not take the fallback, which would hand back that same threshold.
        /// </summary>
        private static string SpectrumTitle(IEnumerable<(string Accession, string Value)> cvParams)
        {
            if (cvParams == null)
            {
                return null;
            }

            var all = cvParams.ToList();
            int title = all.FindIndex(cv => cv.Accession == SpectrumTitleAccession || cv.Accession == ObsoleteSpectrumTitleAccession);
            if (title < 0)
            {
                return all.Select(cv => cv.Value).FirstOrDefault();
            }

            return string.IsNullOrEmpty(all[title].Value) ? null : all[title].Value;
        }

        public string Ms2SpectrumID(int sirIndex)
        {
            string ms2id = null;
            if (dd110 != null)
            {
                var format = dd110.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam;
                var result = dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex];
                if (IsFileFormat(format.accession, format.name, ThermoRawFormatAccession, "Thermo RAW format")
                    || IsFileFormat(format.accession, format.name, MzmlFormatAccession, "mzML format"))
                {
                    ms2id = result.spectrumID;
                }
                else if (IsFileFormat(format.accession, format.name, MascotMgfFormatAccession, "Mascot MGF format"))
                {
                    ms2id = SpectrumTitle(result.cvParam?.Select(cv => (cv.accession, cv.value)));
                }
            }
            else if (dd111 != null)
            {
                var format = dd111.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam;
                var result = dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex];
                if (IsFileFormat(format.accession, format.name, ThermoRawFormatAccession, "Thermo RAW format")
                    || IsFileFormat(format.accession, format.name, MzmlFormatAccession, "mzML format"))
                {
                    ms2id = result.spectrumID;
                }
                else if (IsFileFormat(format.accession, format.name, MascotMgfFormatAccession, "Mascot MGF format"))
                {
                    ms2id = SpectrumTitle(result.cvParam?.Select(cv => (cv.accession, cv.value)));
                }
            }
            else if (dd120 != null)
            {
                var format = dd120.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam;
                var result = dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex];
                if (IsFileFormat(format.accession, format.name, ThermoRawFormatAccession, "Thermo RAW format")
                    || IsFileFormat(format.accession, format.name, MzmlFormatAccession, "mzML format"))
                {
                    ms2id = result.spectrumID;
                }
                else if (IsFileFormat(format.accession, format.name, MascotMgfFormatAccession, "Mascot MGF format"))
                {
                    ms2id = SpectrumTitle(result.cvParam?.Select(cv => (cv.accession, cv.value)));
                }
            }
            else
            {
                var format = dd130.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam;
                var result = dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex];
                if (IsFileFormat(format.accession, format.name, ThermoRawFormatAccession, "Thermo RAW format")
                    || IsFileFormat(format.accession, format.name, MzmlFormatAccession, "mzML format"))
                {
                    ms2id = result.spectrumID;
                }
                else if (IsFileFormat(format.accession, format.name, MascotMgfFormatAccession, "Mascot MGF format"))
                {
                    ms2id = SpectrumTitle(result.cvParam?.Select(cv => (cv.accession, cv.value)));
                }
            }
            return ms2id;
        }

        public float[] MatchedIons(int sirIndex, int siiIndex, int i)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values;
            }
        }

        public int MatchedIonCounts(int sirIndex, int siiIndex, int i)
        {
            if (dd110 != null)
            {
                return dd110.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
            }
            else if (dd111 != null)
            {
                return dd111.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
            }
            else if (dd120 != null)
            {
                return dd120.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
            }
            else
            {
                return dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].Fragmentation[i].FragmentArray[0].values.Length;
            }
        }

        public string ProteinAccession(int sirIndex, int siiIndex)
        {
            string s = null;

            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                string peptideEvidenceRef = 
                    dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef[0].peptideEvidence_ref;
                foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                {
                    if (ok.id.Equals(peptideEvidenceRef))
                    {
                        foreach (var ok2 in dd130.SequenceCollection.DBSequence)
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

        public string ProteinFullName(int sirIndex, int siiIndex)
        {
            string s = "";

            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                foreach (mzIdentML130.Generated.PeptideEvidenceRefType pe 
                    in dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
                    {
                        if (ok.id.Equals(peptideEvidenceRef))
                        {
                            foreach (var ok2 in dd130.SequenceCollection.DBSequence)
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
            return s;
        }

        public string StartResidueInProtein(int sirIndex, int siiIndex)
        {
            string startResidue = "";
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                foreach (mzIdentML130.Generated.PeptideEvidenceRefType pe 
                    in dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
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
            return startResidue;
        }

        public string EndResidueInProtein(int sirIndex, int siiIndex)
        {
            string endResidue = "";
            if (dd110 != null)
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
            else if (dd111 != null)
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
            else if (dd120 != null)
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
            else
            {
                foreach (mzIdentML130.Generated.PeptideEvidenceRefType pe 
                    in dd130.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[sirIndex].SpectrumIdentificationItem[siiIndex].PeptideEvidenceRef)
                {
                    string peptideEvidenceRef = pe.peptideEvidence_ref;
                    foreach (var ok in dd130.SequenceCollection.PeptideEvidence)
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
            return endResidue;
        }

        #endregion Public Methods
    }
}