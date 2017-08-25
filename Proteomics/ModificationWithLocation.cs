using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationWithLocation : Modification
    {
        #region Public Fields

        public static readonly Dictionary<string, TerminusLocalization> terminusLocalizationTypeCodes;
        public readonly IDictionary<string, IList<string>> linksToOtherDbs;
        public readonly List<string> keywords;
        public readonly TerminusLocalization terminusLocalization;
        public readonly ModificationMotif motif;

        #endregion Public Fields

        #region Public Constructors

        static ModificationWithLocation()
        {
            terminusLocalizationTypeCodes = new Dictionary<string, TerminusLocalization>
            {
                { "N-terminal.", TerminusLocalization.NProt }, // Implies protein only, not peptide
                { "C-terminal.", TerminusLocalization.ProtC },
                { "Peptide N-terminal.", TerminusLocalization.NPep },
                { "Peptide C-terminal.", TerminusLocalization.PepC },
                { "Anywhere.", TerminusLocalization.Any },
                { "Protein core.", TerminusLocalization.Any }
            };
        }

        public ModificationWithLocation(string id, string modificationType, ModificationMotif motif, TerminusLocalization terminusLocalization, IDictionary<string, IList<string>> linksToOtherDbs = null, List<string> keywords = null) : base(id, modificationType)
        {
            this.motif = motif;
            this.terminusLocalization = terminusLocalization;

            // Optional
            this.linksToOtherDbs = linksToOtherDbs ?? new Dictionary<string, IList<string>>();
            this.keywords = keywords ?? new List<string>();
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("PP   " + terminusLocalizationTypeCodes.First(b => b.Value.Equals(terminusLocalization)).Key);
            foreach (var nice in linksToOtherDbs)
                foreach (var db in nice.Value)
                    sb.AppendLine("DR   " + nice.Key + "; " + db);
            sb.Append("TG   " + motif);
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithLocation m = o as ModificationWithLocation;
            return m != null
               && base.Equals(m)
               && motif.Equals(m.motif)
               && terminusLocalization == m.terminusLocalization;
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() ^ terminusLocalization.GetHashCode() ^ motif.GetHashCode();
        }

        #endregion Public Methods
    }
}