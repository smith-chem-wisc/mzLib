using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationWithLocation : Modification
    {

        #region Public Fields

        public static readonly Dictionary<string, ModificationSites> modificationTypeCodes;
        public readonly Tuple<string, string> accession;
        public readonly IDictionary<string, IList<string>> linksToOtherDbs;
        public readonly ModificationSites position;
        public readonly ModificationMotif motif;
        public readonly string modificationType;

        #endregion Public Fields

        #region Public Constructors

        static ModificationWithLocation()
        {
            modificationTypeCodes = new Dictionary<string, ModificationSites>();
            modificationTypeCodes.Add("N-terminal.", ModificationSites.NProt); // Implies protein only, not peptide
            modificationTypeCodes.Add("C-terminal.", ModificationSites.ProtC);
            modificationTypeCodes.Add("Peptide N-terminal.", ModificationSites.NPep);
            modificationTypeCodes.Add("Peptide C-terminal.", ModificationSites.PepC);
            modificationTypeCodes.Add("Anywhere.", ModificationSites.Any);
            modificationTypeCodes.Add("Protein core.", ModificationSites.Any);
        }

        public ModificationWithLocation(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites position, IDictionary<string, IList<string>> linksToOtherDbs, string modificationType) : base(id)
        {
            this.accession = accession;
            this.motif = motif;
            this.position = position;
            this.linksToOtherDbs = linksToOtherDbs;
            this.modificationType = modificationType;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("PP   " + modificationTypeCodes.First(b => b.Value.Equals(position)).Key);
            sb.AppendLine("TG   " + motif.Motif);
            if (linksToOtherDbs != null)
                foreach (var nice in linksToOtherDbs)
                    foreach (var db in nice.Value)
                        sb.AppendLine("DR   " + nice.Key + "; " + db);
            sb.Append("MT   " + modificationType);
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithLocation m = o as ModificationWithLocation;
            return m == null ? false :
                
                base.Equals(m)
                && (this.accession == null && ((ModificationWithLocation)m).accession == null
                || this.accession != null && ((ModificationWithLocation)m).accession != null
                && this.accession.Item1 == ((ModificationWithLocation)m).accession.Item1
                && this.accession.Item2 == ((ModificationWithLocation)m).accession.Item2)

                && (this.motif == null && ((ModificationWithLocation)m).motif == null
                || this.motif != null && ((ModificationWithLocation)m).motif != null
                && this.motif.Motif == ((ModificationWithLocation)m).motif.Motif)

                && (this.linksToOtherDbs == null && ((ModificationWithLocation)m).linksToOtherDbs == null
                || this.linksToOtherDbs != null && ((ModificationWithLocation)m).linksToOtherDbs != null
                && this.linksToOtherDbs.Keys.All(a => ((ModificationWithLocation)m).linksToOtherDbs.Keys.Contains(a))
                && this.linksToOtherDbs.Values.SelectMany(x => x).All(a => ((ModificationWithLocation)m).linksToOtherDbs.Values.SelectMany(x => x).Contains(a))
                && ((ModificationWithLocation)m).linksToOtherDbs.Keys.All(a => this.linksToOtherDbs.Keys.Contains(a))
                && ((ModificationWithLocation)m).linksToOtherDbs.Values.SelectMany(x => x).All(a => this.linksToOtherDbs.Values.SelectMany(x => x).Contains(a)))

                && this.modificationType == ((ModificationWithLocation)m).modificationType
                && this.position == ((ModificationWithLocation)m).position;
        }

        public override int GetHashCode()
        {
            int hash = base.GetHashCode() ^ position.GetHashCode();
            hash = hash ^ (modificationType == null ? 0 : modificationType.GetHashCode());
            hash = hash ^ (accession == null ? 0 : accession.GetHashCode());
            hash = hash ^ (motif == null ? 0 : motif.Motif.GetHashCode());
            if (linksToOtherDbs != null)
            {
                foreach (string a in linksToOtherDbs.Keys) hash = hash ^ (a == null ? 0 : a.GetHashCode());
                foreach (string b in linksToOtherDbs.Values.SelectMany(c => c)) hash = hash ^ (b == null ? 0 : b.GetHashCode());
            }
            return hash;
        }

        #endregion Public Methods

    }
}