﻿using System.Text.RegularExpressions;

namespace Omics.Modifications
{
    public class ModificationMotif
    {
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
    }
}