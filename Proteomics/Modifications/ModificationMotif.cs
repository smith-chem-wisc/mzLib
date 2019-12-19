using System.Linq;
using System.Text.RegularExpressions;

namespace Proteomics
{
    public class ModificationMotif
    {
        private static readonly Regex ModificationMotifRegex = new Regex(@"^[A-Za-z]+$", RegexOptions.Compiled);
        public readonly string MotifString;
        public readonly int FirstIndexOfUpper;

        private ModificationMotif(string motif)
        {
            MotifString = motif;
            FirstIndexOfUpper = MotifString.IndexOf(motif.ToString().First(b => char.IsUpper(b)));
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

        public override bool Equals(object o)
        {
            ModificationMotif m = o as ModificationMotif;
            return m != null
                && m.MotifString == MotifString;
        }

        public override int GetHashCode()
        {
            return MotifString.GetHashCode();
        }

        public override string ToString()
        {
            return MotifString;
        }
    }
}