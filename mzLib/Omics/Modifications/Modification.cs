using Chemistry;
using MassSpectrometry;
using System.Text;

namespace Omics.Modifications
{
    /// <summary>
    /// Represents a modification
    /// Mods.txt format was taken from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt
    /// </summary>
    /// <remarks>
    /// IDENTITY IS MUTABLE BY A SUBCLASS, and anything derived from it must survive that.
    /// <see cref="IdWithMotif"/>, <see cref="OriginalId"/>, <see cref="Target"/> and
    /// <see cref="LocationRestriction"/> have protected setters, and MetaMorpheus's
    /// <c>Glycan : Modification</c> really does assign three of them in its own constructor, AFTER
    /// the base constructor has run. So a value memoised from them on first read can be stale by
    /// second read; <see cref="BlocksCleavage"/> handles that by caching the INPUTS next to the
    /// answer and recomputing when they change, rather than by assuming they cannot. Instances are
    /// shared across threads (the static lists on <see cref="Mods"/>), so a stale answer would be
    /// globally wrong rather than locally wrong -- do not replace that with a plain memo.
    /// </remarks>
    public class Modification : IComparable<Modification>
    {
        public string IdWithMotif { get; protected set; }

        /// <summary>
        /// The name of the Mod. This is what shows up in the full sequence. 
        /// </summary>
        public string OriginalId { get; protected set; }
        public string Accession { get; protected set; }

        /// <summary>
        /// The group the modification belongs to. Determines grouping in MetaMorpheuse drop down selections. (Common Biological, Common Fixed)
        /// </summary>
        public string ModificationType { get; protected set; }
        public string FeatureType { get; protected set; }
        public ModificationMotif Target { get; protected set; }

        /// <summary>
        /// Determines where a mod can be placed in an IBioPolymerWithSetMods during digestion. Fixed terminology is stored as strings and found at ModLocationOnPeptideOrProtein and is used throughout the codebase, bit of a mess.  
        /// </summary>
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

        /// <summary>
        /// True when this modification sits on the side chain of a trypsin-family cleavage residue
        /// (Lys/Arg) and neutralises or masks its charge enough that the protease would not cleave
        /// after it -- N6-succinyllysine, N6-acetyllysine and the other epsilon-amine acylations.
        ///
        /// This is a property of the modification alone, and it is only the first of three questions.
        /// Whether the configured PROTEASE cleaves after that residue at all is the second -- see
        /// <see cref="CleavageBlockingModifications.BlocksCleavageBy"/>, which is what digestion
        /// actually consults. Whether the modification invalidates a given peptidoform is the third, a
        /// question of POSITION that digestion decides: an acylated residue that is the protein's own
        /// C-terminus ends a perfectly real peptide, since no cleavage happens there.
        ///
        /// Curated classification; see <see cref="CleavageBlockingModifications"/> for what is in the
        /// set and why the methyl series is excluded. Digestion consults this only when
        /// DigestionParams.RespectCleavageBlockingModifications is set.
        /// </summary>
        public bool BlocksCleavage
        {
            get
            {
                // Classifying costs a lower-casing allocation and a scan of the acyl stems, and the
                // digestion path asks per modification per peptidoform -- millions of times across a
                // search -- so the answer is memoised.
                //
                // The memo cannot assume its inputs are fixed. OriginalId, IdWithMotif, Target and
                // LocationRestriction have protected setters, and MetaMorpheus's Glycan subclass
                // assigns three of them in its own constructor after base construction -- so a plain
                // "compute once" cache can be built from values that no longer hold, and these
                // instances are shared across threads. The cached inputs are therefore stored beside
                // the cached answer and compared by reference on every read: three reference
                // comparisons, far cheaper than the scan, and correct under mutation.
                //
                // The snapshot is a single immutable object behind one reference field, so publishing
                // it is atomic and a read can never see a half-updated pair. The race is benign: two
                // threads racing here compute the same answer from the same inputs.
                string id = OriginalId ?? IdWithMotif;
                ModificationMotif target = Target;
                string locationRestriction = LocationRestriction;

                BlocksCleavageClassification cached = _blocksCleavage;
                if (cached is null
                    || !ReferenceEquals(cached.Id, id)
                    || !ReferenceEquals(cached.Target, target)
                    || !ReferenceEquals(cached.LocationRestriction, locationRestriction))
                {
                    cached = new BlocksCleavageClassification(id, target, locationRestriction,
                        CleavageBlockingModifications.NeutralizesCleavageResidue(this));
                    _blocksCleavage = cached;
                }

                return cached.Answer;
            }
        }

        /// <summary>
        /// A memoised <see cref="BlocksCleavage"/> answer together with the three inputs it was
        /// computed from, so the answer can be invalidated when a subclass changes them.
        /// </summary>
        private sealed class BlocksCleavageClassification
        {
            internal BlocksCleavageClassification(string id, ModificationMotif target, string locationRestriction, bool answer)
            {
                Id = id;
                Target = target;
                LocationRestriction = locationRestriction;
                Answer = answer;
            }

            internal readonly string Id;
            internal readonly ModificationMotif Target;
            internal readonly string LocationRestriction;
            internal readonly bool Answer;
        }

        // Not serialized anywhere: MetaMorpheus's NetSerializer peptide index is field-based and
        // schema-rigid, but PeptideWithSetModifications keeps its modification dictionary
        // [NonSerialized] and rebuilds it from the full sequence, so no Modification reaches that
        // cache and its on-disk layout is unaffected by this field.
        private BlocksCleavageClassification _blocksCleavage;

        public virtual bool ValidModification
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
                case "Protein core.":
                case "Protein core":
                    return "Anywhere.";

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
            int cf = ChemicalFormula?.GetHashCode() ?? 1;
            return id.GetHashCode() ^ mt.GetHashCode() ^ cf;
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
