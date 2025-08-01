﻿using Chemistry;
using MassSpectrometry;
using System.Text;

namespace Omics.Modifications
{
    /// <summary>
    /// Represents a modification
    /// Mods.txt format was taken from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt
    /// </summary>
    public class Modification : IComparable<Modification>
    {
        public string IdWithMotif { get; protected set; }
        public string OriginalId { get; protected set; }
        public string Accession { get; protected set; }
        public string ModificationType { get; protected set; }
        public string FeatureType { get; protected set; }
        public ModificationMotif Target { get; protected set; }
        public string LocationRestriction { get; protected set; }
        public ChemicalFormula ChemicalFormula { get; protected set; }
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

        public Dictionary<string, IList<string>> DatabaseReference { get; protected set; }
        public Dictionary<string, IList<string>> TaxonomicRange { get; protected set; }
        public List<string> Keywords { get; protected set; }
        public Dictionary<DissociationType, List<double>> NeutralLosses { get; protected set; }
        public Dictionary<DissociationType, List<double>> DiagnosticIons { get; protected set; }
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
            if (_originalId != null)
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
                else if (_target != null)
                {
                    this.IdWithMotif = _originalId + " on " + _target.ToString();
                    this.OriginalId = _originalId;
                }
                else
                {
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
                case "C-terminal.":
                case "Peptide N-terminal.":
                case "Peptide C-terminal.":
                case "Anywhere.":
                case "3'-terminal.":
                case "5'-terminal.":
                case "Oligo 3'-terminal.":
                case "Oligo 5'-terminal.":
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
            if (this.NeutralLosses != null)
            {
                if (this.NeutralLosses.Count != 0)
                {
                    List<DissociationType> allDissociationTypes = this.NeutralLosses.Keys.ToList();
                    allDissociationTypes.Sort();

                    foreach (DissociationType dissociationType in allDissociationTypes)
                    {
                        StringBuilder myLine = new StringBuilder();
                        myLine.Append("NL   ");

                        List<double> myValues = new List<double>(this.NeutralLosses[dissociationType]);
                        myValues.Sort();
                        for (int i = 0; i < myValues.Count; i++)
                        {
                            myLine.Append(dissociationType + ":" + ClassExtensions.RoundedDouble(myValues[i]));
                            if (i < myValues.Count - 1)
                                myLine.Append(" or ");
                        }

                        sb.AppendLine(myLine.ToString());
                    }
                }
            }
            if (this.DiagnosticIons != null)
            {
                if (this.DiagnosticIons.Count != 0)
                {
                    List<DissociationType> allDissociationTypes = this.DiagnosticIons.Keys.ToList();
                    allDissociationTypes.Sort();

                    foreach (DissociationType dissociationType in allDissociationTypes)
                    {
                        StringBuilder myLine = new StringBuilder();
                        myLine.Append("DI   ");

                        List<double> myValues = new List<double>(this.DiagnosticIons[dissociationType]);
                        myValues.Sort();
                        for (int i = 0; i < myValues.Count; i++)
                        {
                            myLine.Append(dissociationType + ":" + ClassExtensions.RoundedDouble(myValues[i]));
                            if (i < myValues.Count - 1)
                                myLine.Append(" or ");
                        }

                        sb.AppendLine(myLine.ToString());
                    }
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


        // Used in the sorted sets for variable mod generation to ensure that modifications are consistently ordered
        // UniProt annotations also contain an evidence level. Future work could include this in the ordering of modifications for digestion. 
        public int CompareTo(Modification? other)
        {
            if (other == null) return 1;

            int idComparison = string.Compare(this.IdWithMotif, other.IdWithMotif, StringComparison.Ordinal);
            if (idComparison != 0) return idComparison;

            int typeComparison = string.Compare(this.ModificationType, other.ModificationType, StringComparison.Ordinal);
            if (typeComparison != 0) return typeComparison;

            int locRestrictionComparison = string.Compare(this.LocationRestriction, other.LocationRestriction, StringComparison.Ordinal);
            if (locRestrictionComparison != 0) return locRestrictionComparison;

            return Nullable.Compare(this.MonoisotopicMass, other.MonoisotopicMass);
        }
    }
}