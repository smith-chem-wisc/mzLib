using System.Text.RegularExpressions;

namespace Omics.Modifications
{
    public class ModificationMotif
    {
        private static Dictionary<string, ModificationMotif> MotifCache = []; // Cache for ModificationMotifs to reduce memory overhead by using shared references. 
        private static readonly Regex ModificationMotifRegex = new Regex(@"^[A-Za-z]+$", RegexOptions.Compiled);
        private readonly string motifString;

        private ModificationMotif(string motif)
        {
            motifString = motif;
        }

        /// <summary>
        /// Only upper and lower case letters allowed, must have a single upper case letter
        /// </summary>
        /// <param name="motifString"></param>
        /// <param name="motif"></param>
        /// <returns></returns>
        public static bool TryGetMotif(string motifString, out ModificationMotif motif)
        {
            motif = null!;
            if (ModificationMotifRegex.IsMatch(motifString) && motifString.Count(b => char.IsUpper(b)) == 1)
            {
                // Check Motif Cache to see if we have generated this one before. 
                if (MotifCache.TryGetValue(motifString, out motif))
                {
                    return true;
                }
                else
                {
                    motif = new ModificationMotif(motifString);
                    MotifCache[motifString] = motif;
                    return true;
                }
            }
            return false;
        }

        public override string ToString()
        {
            return motifString;
        }
    }
}