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

        public override bool Equals(Modification m)
        {
            return base.Equals(m)
                && m as ModificationWithLocation != null
                
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

        public override int GetCustomHashCode()
        {
            int hash = base.GetCustomHashCode() + sum_string_chars(position.ToString()) + sum_string_chars(modificationType);
            if (accession != null) hash += sum_string_chars(accession.Item1) + sum_string_chars(accession.Item2);
            if (motif != null) hash += sum_string_chars(motif.Motif);
            if (linksToOtherDbs != null)
            {
                foreach (string a in linksToOtherDbs.Keys) hash += sum_string_chars(a);
                foreach (string b in linksToOtherDbs.Values.SelectMany(c => c)) hash += sum_string_chars(b);
            }
            return hash;
        }

        #endregion Public Methods

        #region Private Methods

        private int sum_string_chars(string s)
        {
            if (s == null) return 0;
            int sum = int.MinValue + s.ToCharArray().Sum(c => 37 * c);
            return sum != 0 ? sum : sum + 1;
        }

        #endregion Private Methods

    }
}