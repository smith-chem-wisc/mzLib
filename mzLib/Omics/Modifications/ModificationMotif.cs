using System.Text.RegularExpressions;

namespace Omics.Modifications
{
    public class ModificationMotif
    {
        private static readonly Regex ModificationMotifRegex = new Regex(@"^[A-Za-z]+$", RegexOptions.Compiled);
        private readonly string motifString;
        public Regex ModficationMotifPattern;

        private ModificationMotif(string motif)
        {
            motifString = motif;
            string modifiedSequencePattern = Regex.Escape(motif) + @"\[[^\]]+\]";
            ModficationMotifPattern = new Regex(modifiedSequencePattern);
        }

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
        
        public override string ToString()
        {
            return motifString;
        }
        public override int GetHashCode()
        {
            return motifString.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (obj is ModificationMotif other)
            {
                return string.Equals(motifString, other.motifString, StringComparison.Ordinal);
            }
            return false;
        }
    }
}