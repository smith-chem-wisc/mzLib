using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationGeneral
    {
        public string Id { get; private set; }
        public string Accession { get; private set; }
        public string ModificationType { get; private set; }
        public string FeatureType { get; private set; }
        public ModificationMotif Target { get; private set; }
        public string Position { get; private set; }
        public ChemicalFormula ChemicalFormula { get; private set; }
        private double? monoisotopicMass = null;

        public double? MonoisotopicMass
        {
            get
            {
                return RoundedDouble(monoisotopicMass);
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

        public bool ValidModification
        {
            get { return (this.Id != null && (this.ChemicalFormula != null || this.MonoisotopicMass != null) && this.Position != "Unassigned." && this.ModificationType != null && this.FeatureType != "CROSSLINK"); }
            private set { ValidModification = value; }
        }

        public ModificationGeneral(string _Id = null, string _Accession = null, string _ModificationType = null, string _FeatureType = null, ModificationMotif _Target = null, string _Position = "Unassigned.", ChemicalFormula _ChemicalFormula = null, double? _MonoisotopicMass = null, Dictionary<string, IList<string>> _DatabaseReference = null, Dictionary<string, IList<string>> _TaxonomicRange = null, List<string> _Keywords = null, Dictionary<DissociationType, List<double>> _NeutralLosses = null, Dictionary<DissociationType, List<double>> _DiagnosticIons = null, string _FileOrigin = null)
        {
            this.Id = _Id;
            this.Accession = _Accession;
            this.ModificationType = _ModificationType;
            this.FeatureType = _FeatureType;
            this.Target = _Target;
            this.Position = ModLocationOnPeptideOrProtein(_Position);
            this.ChemicalFormula = _ChemicalFormula;
            this.MonoisotopicMass = _MonoisotopicMass;
            this.DatabaseReference = _DatabaseReference;
            this.TaxonomicRange = _TaxonomicRange;
            this.Keywords = _Keywords;
            this.NeutralLosses = _NeutralLosses;
            this.DiagnosticIons = _DiagnosticIons;
            this.FileOrigin = _FileOrigin;
        }

        public static string ModLocationOnPeptideOrProtein(string modLocation)
        {
            List<string> locationList = new List<string>
            {
                "N-terminal.",
                "C-terminal.",
                "Peptide N-terminal.",
                "Peptide C-terminal.",
                "Anywhere."
            };
            if (locationList.Contains(modLocation))
            {
                return modLocation;
            }
            else
            {
                return "Unassigned.";
            }
        }

        public double? RoundedDouble(double? myNumber)
        {
            if (myNumber != null)
            {
                myNumber = Math.Round((double)myNumber, 9, MidpointRounding.AwayFromZero);
            }
            return myNumber;
        }

        public override bool Equals(object o)
        {
            ModificationGeneral m = o as ModificationGeneral;
            return o != null
                && m.ToString() == this.ToString();
        }

        public override int GetHashCode()
        {
            return this.ToString().GetHashCode();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            if (this.Id != null)
            { sb.AppendLine("ID   " + this.Id); }
            if (this.Accession != null)
            { sb.AppendLine("AC   " + this.Accession); }
            if (this.ModificationType != null)
            { sb.AppendLine("MT   " + this.ModificationType); }
            if (this.FeatureType != null)
            { sb.AppendLine("FT   " + this.FeatureType); }
            if (this.Target != null)
            { sb.AppendLine("TG   " + this.Target); } // at this stage, each mod has only one target though many may have the same Id
            if (this.Position != null)
            { sb.AppendLine("PP   " + this.Position); }
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
                            myLine.Append(myKey + ":" + RoundedDouble(myValues[i]));
                            if (i < myValues.Count)
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
                            myLine.Append(myKey + ":" + RoundedDouble(myValues[i]));
                            if (i < myValues.Count)
                                myLine.Append(" or ");
                        }
                    }
                    sb.Append(myLine);
                }
            }

            if (this.Keywords != null)
            {
                if (this.Keywords.Count != 0)
                { sb.AppendLine("KW   " + String.Join(" or ", this.Keywords.ToList().OrderBy(b => b))); }
            }

            return sb.ToString();
        }

        public string ModificationErrorsToString() //reports errors in required fields.
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(this.ToString());

            if (this.Id == null)
                sb.AppendLine("#Required field ID missing or malformed. Current value = " + this.Id);
            if (this.ModificationType == null)
                sb.AppendLine("#Required field MT missing or malformed. Current value = " + this.ModificationType);
            if (this.Position == null)
                sb.AppendLine("#Required field PP missing or malformed. Current value = " + this.Position + ".");
            if (this.ChemicalFormula == null && this.MonoisotopicMass == null)
                sb.AppendLine("#Required fields CF and MM are both missing or malformed. One of those two fields must be provided.");
            sb.Append("#This modification can be found in file " + this.FileOrigin);

            return sb.ToString();
        }
    }
}