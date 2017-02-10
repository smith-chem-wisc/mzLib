using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationWithLocation : Modification
    {

        #region Public Fields

        public readonly Tuple<string, string> accession;
        public readonly IDictionary<string, IList<string>> linksToOtherDbs;
        public readonly ModificationSites position;
        public readonly ModificationMotif motif;
        public readonly string database;

        #endregion Public Fields

        #region Private Fields

        public static readonly Dictionary<string, ModificationSites> modificationTypeCodes;

        #endregion Private Fields

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

        public ModificationWithLocation(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites position, IDictionary<string, IList<string>> linksToOtherDbs, string database) : base(id)
        {
            this.accession = accession;
            this.motif = motif;
            this.position = position;
            this.linksToOtherDbs = linksToOtherDbs;
            this.database = database;
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
            sb.Append("DB   " + database);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}