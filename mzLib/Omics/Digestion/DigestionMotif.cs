using System.Text.RegularExpressions;
using MzLibUtil;

namespace Omics.Digestion
{
    public class DigestionMotif
    {
        private static char[] B = new char[] { 'D', 'N' };
        private static char[] J = new char[] { 'I', 'L' };
        private static char[] Z = new char[] { 'E', 'Q' };

        public readonly string InducingCleavage;
        public readonly string PreventingCleavage;
        public readonly int CutIndex;
        public readonly string ExcludeFromWildcard;

        public DigestionMotif(string inducingCleavage, string preventingCleavage, int cutIndex, string excludeFromWildcard)
        {
            this.InducingCleavage = inducingCleavage;
            this.PreventingCleavage = preventingCleavage;
            this.CutIndex = cutIndex;
            this.ExcludeFromWildcard = excludeFromWildcard;
        }

        // parsing cleavage rules syntax
        public static List<DigestionMotif> ParseDigestionMotifsFromString(string motifsString)
        {
            motifsString = motifsString.Replace("\"", string.Empty).Replace(" ", string.Empty);

            // throws exception if non-supported characters are used
            if (Regex.Match(motifsString, @"[^a-zA-Z0-9|,[\]{}]+").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. The digestion motif can only contain letters and {}[]|");
            }
            // throws exception if user attempts separate multiple preventing cleavages using commas
            if (Regex.Match(motifsString, @"\[([\w]*,+[\w]*)*\]").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. Please create a separate motif for each sequence preventing cleavage (comma separated).");
            }
            // throws exception if user attempts separate multiple wildcard exclusions
            if (Regex.Match(motifsString, @"\{([\w]*,+[\w]*)*\}").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. Please create a separate motif for each wildcard exclusion (comma separated).");
            }

            string[] motifStrings = motifsString.Split(',');
            var motifs = new List<DigestionMotif>();

            for (int i = 0; i < motifStrings.Length; i++)
            {
                string motifString = motifStrings[i];
                motifs.Add(ParseDigestionMotifFromString(motifString));
            }
            return motifs;
        }

        private static DigestionMotif ParseDigestionMotifFromString(string motifString)
        {
            string inducingCleavage;
            string preventingCleavage = null;
            string excludingWC = null;
            int cutIndex = 0;

            if (motifString.Contains("{") && !motifString.Contains("}")
                || !motifString.Contains("{") && motifString.Contains("}")
                || motifString.Contains("[") && !motifString.Contains("]")
                || !motifString.Contains("[") && motifString.Contains("]"))
            {
                throw new MzLibException("Unrecognized protease syntax. Please close any brackets used.");
            }

            // find preventing cleavage
            if (motifString.Contains("["))
            {
                int start = motifString.IndexOf("[") + 1;
                int end = motifString.IndexOf("]");

                preventingCleavage = motifString.Substring(start, end - start);
                motifString = Regex.Replace(motifString, @"\[[a-zA-Z]+\]", string.Empty);
            }

            // finds wildcard exceptions
            if (motifString.Contains("{"))
            {
                int start = motifString.IndexOf("{") + 1;
                int end = motifString.IndexOf("}");

                excludingWC = motifString.Substring(start, end - start);
                if (Regex.Matches(motifString.ToUpper(), "X").Count != excludingWC.Length)
                {
                    throw new MzLibException("Unrecognized protease syntax. Please have equal number of wildcards for multi-letter wildcard exclusions.");
                }
                motifString = Regex.Replace(motifString, @"\{[a-zA-Z]+\}", string.Empty);
            }

            // finds motif cut index
            for (int j = 0; j < motifString.Length; j++)
            {
                if (motifString[j] == '|')
                {
                    cutIndex = j;
                    break;
                }
            }

            motifString = motifString.Replace("|", string.Empty);
            inducingCleavage = motifString;

            return new DigestionMotif(inducingCleavage, preventingCleavage, cutIndex, excludingWC);
        }

        public (bool, bool) Fits(string sequence, int location)
        {
            bool fits = true;
            char currentResidue;
            int m;

            // check for inducing cleavage
            for (m = 0; m < InducingCleavage.Length && fits; m++) // handle patterns
            {
                if (location + m >= sequence.Length)
                {
                    fits = false;
                }
                else
                {
                    currentResidue = sequence[location + m];
                    if (!MotifMatches(InducingCleavage[m], currentResidue))
                    {
                        fits = false;
                    }
                }
            }

            bool prevents = false;
            // check for preventing cleavage
            if (fits && PreventingCleavage != null)
            {
                prevents = true;
                for (int n = 0; n < PreventingCleavage.Length && prevents; n++)
                {
                    if (location + m + n >= sequence.Length || location - PreventingCleavage.Length + 1 + n < 0)
                    {
                        prevents = false;
                    }
                    else
                    {
                        currentResidue = CutIndex != 0 ? sequence[location + m + n] : sequence[location - PreventingCleavage.Length + 1 + n];
                        if (!MotifMatches(PreventingCleavage[n], currentResidue))
                        {
                            prevents = false;
                        }
                    }
                }

                fits = prevents ? false : true;
            }

            return (fits, prevents);
        }

        /// <summary>
        /// True when this motif severs the bond C-TERMINAL to <paramref name="residue"/> -- that is,
        /// when the residue whose own C-side bond the cut index falls after matches. Trypsin's "K|"
        /// and "R|" report true for K and R; Asp-N's "|D" and Lys-N's "|K" report false for every
        /// residue, because they cut N-terminal to their recognition residue and sever nothing after it.
        /// </summary>
        /// <remarks>
        /// Residue-level only: a preventing-cleavage rule (trypsin|P's "K[P]|") is deliberately not
        /// consulted, because that rule depends on the sequence context of a particular site and this
        /// question is about the protease alone. Ambiguity codes and the wildcard are honoured through
        /// the same matcher digestion itself uses, so "X|" reports true for every residue.
        ///
        /// This is exactly the P1 subsite question, so it delegates to
        /// <see cref="NonPrimeSubsiteAccepts"/> rather than deriving InducingCleavage[CutIndex - 1] a
        /// second time. The two are equivalent by construction: a CutIndex below 1 or above the
        /// recognition sequence's length puts P1 outside the motif, which is the same rejection the
        /// explicit bounds check used to make.
        /// </remarks>
        public bool CleavesCTerminalTo(char residue) => NonPrimeSubsiteAccepts(1, residue);

        /// <summary>
        /// Where this motif's recognition sequence sits relative to the bond it severs. See
        /// <see cref="CleavageSide"/> for why this has three values and not two -- a motif may straddle
        /// the bond, naming residues on both sides, and four such motifs ship today.
        /// </summary>
        public CleavageSide Side =>
            CutIndex <= 0 ? CleavageSide.NTerminal
            : CutIndex >= (InducingCleavage?.Length ?? 0) ? CleavageSide.CTerminal
            : CleavageSide.Straddling;

        /// <summary>
        /// The residue this motif requires at the non-prime subsite P<paramref name="position"/> --
        /// P1 being the residue immediately BEFORE the cut, P2 the one before that, in the
        /// Schechter-Berger convention. The null character when the motif says nothing there, either
        /// because the subsite lies outside the recognition sequence or because
        /// <paramref name="position"/> is not a subsite number.
        /// </summary>
        /// <remarks>
        /// Returns the motif character as written, so it may be an ambiguity code or the wildcard 'X'
        /// rather than a literal residue -- StcE-trypsin's "TX|T" answers 'X' at P1 and 'T' at P2. Use
        /// <see cref="NonPrimeSubsiteAccepts"/> to ask whether a given residue satisfies the subsite,
        /// which honours those codes; comparing this char directly does not.
        /// </remarks>
        public char NonPrimeSubsite(int position) =>
            position < 1 ? '\0' : SubsiteAt(CutIndex - position);

        /// <summary>
        /// The residue this motif requires at the prime subsite P<paramref name="position"/>' --
        /// P1' being the residue immediately AFTER the cut. The null character when the motif says
        /// nothing there. Asp-N's "|D" answers 'D' at P1'.
        /// </summary>
        /// <remarks>
        /// As for <see cref="NonPrimeSubsite"/>, the character is the motif's own, ambiguity codes
        /// included. Prefer <see cref="PrimeSubsiteAccepts"/> for matching.
        /// </remarks>
        public char PrimeSubsite(int position) =>
            position < 1 ? '\0' : SubsiteAt(CutIndex + position - 1);

        /// <summary>
        /// True when <paramref name="residue"/> satisfies this motif's requirement at the non-prime
        /// subsite P<paramref name="position"/>. False when the motif constrains nothing there -- an
        /// unconstrained subsite is not a match, because the caller is asking whether the motif DEMANDS
        /// this residue, not whether it tolerates it.
        /// </summary>
        public bool NonPrimeSubsiteAccepts(int position, char residue) =>
            SubsiteAccepts(NonPrimeSubsite(position), residue);

        /// <summary>
        /// True when <paramref name="residue"/> satisfies this motif's requirement at the prime subsite
        /// P<paramref name="position"/>'. False when the motif constrains nothing there.
        /// </summary>
        public bool PrimeSubsiteAccepts(int position, char residue) =>
            SubsiteAccepts(PrimeSubsite(position), residue);

        /// <summary>
        /// The motif character at a zero-based index into the recognition sequence, or the null
        /// character when the index falls outside it.
        /// </summary>
        private char SubsiteAt(int index) =>
            InducingCleavage is null || index < 0 || index >= InducingCleavage.Length
                ? '\0'
                : InducingCleavage[index];

        /// <summary>
        /// Matches a subsite's motif character against a sequence residue through the same matcher
        /// digestion itself uses, so the wildcard and the B/J/Z ambiguity codes are honoured. The null
        /// character -- an unconstrained subsite -- never matches, and must be rejected here rather
        /// than handed to <see cref="MotifMatches"/>, which would compare it literally.
        /// </summary>
        private bool SubsiteAccepts(char motifResidue, char residue) =>
            motifResidue != '\0' && MotifMatches(motifResidue, residue);

        private bool MotifMatches(char motifChar, char sequenceChar)
        {
            return motifChar.Equals('X') && !sequenceChar.ToString().Equals(ExcludeFromWildcard)
                || motifChar.Equals(sequenceChar)
                || motifChar.Equals('B') && B.Contains(sequenceChar)
                || motifChar.Equals('J') && J.Contains(sequenceChar)
                || motifChar.Equals('Z') && Z.Contains(sequenceChar);
        }
    }
}