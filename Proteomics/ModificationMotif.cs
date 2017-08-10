using System.Linq;
using System.Text.RegularExpressions;

namespace Proteomics
{
    public class ModificationMotif
    {
        #region Private Fields

        private static readonly Regex ModificationMotifRegex = new Regex(@"^[A-Za-z]+$", RegexOptions.Compiled);

        #endregion Private Fields

        #region Private Constructors

        private ModificationMotif(string motif)
        {
            Motif = motif;
        }

        #endregion Private Constructors

        #region Public Properties

        public string Motif { get; }

        #endregion Public Properties

        #region Public Methods

        /// <summary>
        /// Only upper and lower case letters allowed, must have a single upper case letter
        /// </summary>
        /// <param name="motifString"></param>
        /// <param name="motif"></param>
        /// <returns></returns>
        public static bool TryGetMotif(string motifString, out ModificationMotif motif)
        {
            motif = null;
            if (ModificationMotifRegex.IsMatch(motifString) && motifString.Count(b => char.IsUpper(b)) == 1)
            {
                motif = new ModificationMotif(motifString);
                return true;
            }
            return false;
        }

        #endregion Public Methods
    }
}