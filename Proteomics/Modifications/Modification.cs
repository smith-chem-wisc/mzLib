using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class Modification
    {
        public string IdWithMotif { get; private set; }
        public string OriginalId { get; private set; }
        public string Accession { get; private set; }
        public string ModificationType { get; private set; }
        public string FeatureType { get; private set; }
        public ModificationMotif Target { get; private set; }
        public string LocationRestriction { get; private set; }
        public ChemicalFormula ChemicalFormula { get; private set; }
        private double? monoisotopicMass = null;

        public double? MonoisotopicMass
        {
            get
            {
                return ClassExtensions.RoundedDouble(monoisotopicMass);
            }
            private set
            {
                monoisotopicMass = value;
            }
        }

        public Dictionary<string, IList<string>> DatabaseReference { get; private set; }
        public Dictionary<string, IList<string>> TaxonomicRange { get; private set; }
        public List<string> Keywords { get; private set; }
        public Dictionary<DissociationType, List<double>> NeutralLosses { get; private set; }
        public Dictionary<DissociationType, List<double>> DiagnosticIons { get; private set; }
        public string FileOrigin { get; private set; }
        protected const double tolForEquality = 1e-9;

        public bool ValidModification
        {
            get
            {
                return this.IdWithMotif != null
                       && (this.ChemicalFormula != null || this.MonoisotopicMass != null)
                       && this.Target != null
                       && this.LocationRestriction != "Unassigned."
                       && this.ModificationType != null
                       && this.FeatureType != "CROSSLINK"
                       && !this.ModificationType.Contains(':');
            }
        }

        public Modification(string _originalId = null, string _accession = null, string _modificationType = null, string _featureType = null,
            ModificationMotif _target = null, string _locationRestriction = "Unassigned.", ChemicalFormula _chemicalFormula = null,
            double? _monoisotopicMass = null, Dictionary<string, IList<string>> _databaseReference = null,
            Dictionary<string, IList<string>> _taxonomicRange = null, List<string> _keywords = null,
            Dictionary<DissociationType, List<double>> _neutralLosses = null, Dictionary<DissociationType, List<double>> _diagnosticIons = null,
            string _fileOrigin = null)
        {

            if (_originalId != null && _target != null)
            {
                if (_originalId.Contains(" on "))
                {
                    this.IdWithMotif = _originalId;
                    this.OriginalId = _originalId.Split(new[] { " on " }, StringSplitOptions.None)[0];
                }
                else if (_originalId.Contains(" of "))
                {
                    this.IdWithMotif = _originalId.Replace(" of ", " on ");
                    this.OriginalId = _originalId.Split(new[] { " of ", " on " }, StringSplitOptions.None)[0];

                }
                else
                {
                    this.IdWithMotif = _originalId + " on " + _target.ToString();
                    this.OriginalId = _originalId;
                }
            }

            this.Accession = _accession;
            this.ModificationType = _modificationType;
            this.FeatureType = _featureType;
            this.Target = _target;
            this.LocationRestriction = ModLocationOnPeptideOrProtein(_locationRestriction);
            this.ChemicalFormula = _chemicalFormula;
            this.MonoisotopicMass = _monoisotopicMass;
            this.DatabaseReference = _databaseReference;
            this.TaxonomicRange = _taxonomicRange;
            this.Keywords = _keywords;
            this.NeutralLosses = _neutralLosses;
            this.DiagnosticIons = _diagnosticIons;
            this.FileOrigin = _fileOrigin;

            if (this.MonoisotopicMass == null && this.ChemicalFormula != null)
            {
                this.MonoisotopicMass = this.ChemicalFormula.MonoisotopicMass;
            }
        }

        public static string ModLocationOnPeptideOrProtein(string _locationRestriction)
        {
            switch (_locationRestriction)
            {
                case "N-terminal.":
                    return _locationRestriction;

                case "C-terminal.":
                    return _locationRestriction;

                case "Peptide N-terminal.":
                    return _locationRestriction;

                case "Peptide C-terminal.":
                    return _locationRestriction;

                case "Anywhere.":
                    return _locationRestriction;

                default:
                    return "Unassigned.";
            }
        }

        public override bool Equals(object o)
        {
            Modification m = o as Modification;
            return o != null
                && IdWithMotif == m.IdWithMotif
                && OriginalId == m.OriginalId
                && ModificationType == m.ModificationType
                && (MonoisotopicMass == m.MonoisotopicMass
                    || MonoisotopicMass != null && m.MonoisotopicMass != null && Math.Abs((double)m.MonoisotopicMass - (double)MonoisotopicMass) < tolForEquality);
        }

        public override int GetHashCode()
        {
            string id = IdWithMotif ?? OriginalId ?? string.Empty;
            string mt = ModificationType ?? string.Empty;
            return id.GetHashCode() ^ mt.GetHashCode();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            if (this.IdWithMotif != null)
            { sb.AppendLine("ID   " + this.IdWithMotif); }
            if (this.Accession != null)
            { sb.AppendLine("AC   " + this.Accession); }
            if (this.ModificationType != null)
            { sb.AppendLine("MT   " + this.ModificationType); }
            if (this.FeatureType != null)
            { sb.AppendLine("FT   " + this.FeatureType); }
            if (this.Target != null)
            { sb.AppendLine("TG   " + this.Target); } // at this stage, each mod has only one target though many may have the same Id
            if (this.LocationRestriction != null)
            { sb.AppendLine("PP   " + this.LocationRestriction); }
            if (this.ChemicalFormula != null)
            { sb.AppendLine("CF   " + this.ChemicalFormula.Formula); }
            if (this.MonoisotopicMass != null)
            { sb.AppendLine("MM   " + this.MonoisotopicMass); }
            if (this.DatabaseReference != null)
            {
                if (this.DatabaseReference.Count != 0)
                {
                    List<string> myKeys = new List<string>(this.DatabaseReference.Keys);
                    myKeys.Sort();
                    foreach (string myKey in myKeys)
                    {
                        List<string> myValues = new List<string>(this.DatabaseReference[myKey]);
                        myValues.Sort();
                        foreach (string myValue in myValues)
                        {
                            sb.AppendLine("DR   " + myKey + "; " + myValue);
                        }
                    }
                }
            }
            if (this.TaxonomicRange != null)
            {
                if (this.TaxonomicRange.Count != 0)
                {
                    List<string> myKeys = new List<string>(this.TaxonomicRange.Keys);
                    myKeys.Sort();
                    foreach (string myKey in myKeys)
                    {
                        List<string> myValues = new List<string>(this.TaxonomicRange[myKey]);
                        myValues.Sort();
                        foreach (string myValue in myValues)
                        {
                            sb.AppendLine("TR   " + myKey + "; " + myValue);
                        }
                    }
                }
            }
            if (this.Keywords != null)
            {
                if (this.Keywords.Count != 0)
                {
                    sb.Append("KW   " + String.Join(" or ", this.Keywords));
                }
            }
            if (this.NeutralLosses != null)
            {
                if (this.NeutralLosses.Count != 0)
                {
                    StringBuilder myLine = new StringBuilder();
                    myLine.Append("NL   ");
                    List<DissociationType> myKeys = new List<DissociationType>(this.NeutralLosses.Keys);
                    myKeys.Sort();
                    foreach (DissociationType myKey in myKeys)
                    {
                        List<double> myValues = new List<double>(this.NeutralLosses[myKey]);
                        myValues.Sort();
                        for (int i = 0; i < myValues.Count; i++)
                        {
                            myLine.Append(myKey + ":" + ClassExtensions.RoundedDouble(myValues[i]));
                            if (i < myValues.Count - 1)
                                myLine.Append(" or ");
                        }
                    }
                    sb.Append(myLine);
                }
            }
            if (this.DiagnosticIons != null)
            {
                if (this.DiagnosticIons.Count != 0)
                {
                    StringBuilder myLine = new StringBuilder();
                    myLine.Append("DI   ");
                    List<DissociationType> myKeys = new List<DissociationType>(this.DiagnosticIons.Keys);
                    myKeys.Sort();
                    foreach (DissociationType myKey in myKeys)
                    {
                        List<double> myValues = new List<double>(this.DiagnosticIons[myKey]);
                        myValues.Sort();

                        for (int i = 0; i < myValues.Count; i++)
                        {
                            myLine.Append(myKey + ":" + ClassExtensions.RoundedDouble(myValues[i]));

                            if (i < myValues.Count - 1)
                            {
                                myLine.Append(" or ");
                            }
                        }
                    }
                    sb.Append(myLine);
                }
            }

            if (this.Keywords != null)
            {
                if (this.Keywords.Count != 0)
                {
                    sb.AppendLine("KW   " + String.Join(" or ", this.Keywords.ToList().OrderBy(b => b)));
                }
            }

            return sb.ToString();
        }

        public string ModificationErrorsToString() //reports errors in required fields.
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(this.ToString());

            if (this.IdWithMotif == null)
            {
                sb.AppendLine("#Required field ID missing or malformed. Current value = " + this.IdWithMotif);
            }

            if (this.ModificationType == null)
            {
                sb.AppendLine("#Required field MT missing or malformed. Current value = " + this.ModificationType);
            }

            if (this.LocationRestriction == null)
            {
                sb.AppendLine("#Required field PP missing or malformed. Current value = " + this.LocationRestriction +
                              ".");
            }

            if (this.ChemicalFormula == null && this.MonoisotopicMass == null)
            {
                sb.AppendLine(
                    "#Required fields CF and MM are both missing or malformed. One of those two fields must be provided.");
            }

            if (this.ModificationType != null && this.ModificationType.Contains(':'))
            {
                sb.AppendLine("#Modification type cannot contain ':'!");
            }

            sb.Append("#This modification can be found in file " + this.FileOrigin);

            return sb.ToString();
        }
    }
}