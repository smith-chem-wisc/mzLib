using System;
using MzLibUtil;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Proteomics.ProteolyticDigestion
{
    public class DigestionMotif
    {
        public string InducingCleavage;
        public string PreventingCleavage;
        public int CutIndex;
     
        public DigestionMotif(string inducingCleavage, string preventingCleavage, int cutIndex)
        {
            this.InducingCleavage = inducingCleavage;
            this.PreventingCleavage = preventingCleavage;
            this.CutIndex = cutIndex;
        }

        // parsing cleavage rules syntax
        public static List<DigestionMotif> ParseProteaseFromString(string motifString)
        {
            string inducingCleavage = null;
            string preventingCleavage = null;
            int cutIndex = 0;
            
            string[] sequences = motifString.Replace("\"", "").Replace(" ", string.Empty).Split(',');
            var motifs = new List<DigestionMotif>();

            for (int i = 0; i < sequences.Length; ++i)
            {
                string s = sequences[i];

                // error checking for syntax: missing "|" for cleavage terminus
                if (!String.IsNullOrEmpty(s) && !s.Contains("|"))
                    throw new MzLibException("Please specify cleavage terminus.");

                int j = 0;
                bool done = false;

                // finds protease cut index
                while (!done && j < s.Length)
                {
                    switch (s[j])
                    {
                        case '[':
                            goto case '|';
                        case '|':
                            cutIndex = j;
                            done = true;
                            break;
                    }
                    ++j;
                }
                s = s.Replace("|", "");

                // identify inducing and preventing cleavage sequences 
                if (!(s.Contains("["))) 
                {
                    inducingCleavage = s;
                    preventingCleavage = null;
                }
                else
                {
                    int start = s.IndexOf("[") + 1;
                    int end = s.IndexOf("]");
                    inducingCleavage = s.Substring(0, start - 1);
                    preventingCleavage = s.Substring(start, end - start);
                }

               motifs.Add(new DigestionMotif(inducingCleavage, preventingCleavage, cutIndex));
            }

            return motifs;
        }
    }
}
