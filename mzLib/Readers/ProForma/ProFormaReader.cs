using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// mzLib entry point for parsing a HUPO-PSI ProForma 2.0 proteoform string into a
    /// <see cref="Tdp.ProFormaTerm"/> (Layer 1, lossless string -> term).
    /// Parsing is delegated to the TopDownProteomics reference implementation; this facade
    /// exists so mzLib/MetaMorpheus callers depend on a stable mzLib type, not the SDK directly.
    /// </summary>
    /// <remarks>
    /// A term holds one modification per terminus, and the SDK (0.0.299) neither refuses a second one nor stops
    /// residues after the C-terminal modification: <c>PEPTIDE-[A][B]</c> kept only B, <c>[A][B]-PEPTIDE</c> was misread
    /// as an unlocalized A (after an unlocalized prefix, <c>[P]?[A][B]-PEPTIDE</c>, it read A as unlocalized without an
    /// error), <c>PEPTIDE-[A]K</c> read as PEPTIDEK, and a string with no residue threw an index exception or read as an
    /// empty sequence. Before the SDK runs, this reader refuses those shapes, so it never returns a term that lost or
    /// moved a terminal modification.
    /// </remarks>
    public static class ProFormaReader
    {
        /// <summary>
        /// Parses a ProForma string into its term representation.
        /// </summary>
        /// <param name="proFormaString">A ProForma 2.0 string.</param>
        /// <returns>The parsed <see cref="Tdp.ProFormaTerm"/>.</returns>
        /// <exception cref="ProFormaUnsupportedException">The string carries several modifications on one terminus: valid
        /// ProForma 2.1, but a term cannot hold them.</exception>
        /// <exception cref="Tdp.ProFormaParseException">Thrown when the input is not valid ProForma, including a residue
        /// after the C-terminal modification and a string with no residue.</exception>
        public static Tdp.ProFormaTerm Read(string proFormaString)
        {
            ThrowIfTerminalsCannotBeRead(proFormaString);
            return new Tdp.ProFormaParser().ParseString(proFormaString);
        }

        /// <summary>
        /// Refuses, before the SDK sees them: several modification groups on one terminus; a residue after the C-terminal
        /// modification; a string with a modification but no residue. Unlocalized groups (<c>[Phospho]?</c>,
        /// <c>[Phospho]^2?</c>) are passed over first, so the N-terminal group is the one after the last '?'. Bracket,
        /// brace and angle groups are skipped whole, so brackets inside a modification (a formula with isotopes) are never
        /// read as a terminus. Any shape it does not recognise is left to the SDK.
        /// </summary>
        internal static void ThrowIfTerminalsCannotBeRead(string proForma)
        {
            if (string.IsNullOrEmpty(proForma))
                return;
            int n = proForma.Length, i = 0;
            while (i < n && (proForma[i] == '<' || proForma[i] == '{'))
            {
                i = SkipGroup(proForma, i);
                if (i < 0) return;
            }

            // Groups that a '?' closes are unlocalized; only those after the last '?' can be N-terminal.
            int leading = 0;
            bool modified = false;
            while (i < n)
            {
                char c = proForma[i];
                if (c is '[' or '{')
                {
                    i = SkipGroup(proForma, i);
                    if (i < 0) return;
                    if (c == '[') { leading++; modified = true; }
                }
                else if (c == '^' && leading > 0)
                {
                    do i++; while (i < n && char.IsDigit(proForma[i]));
                }
                else if (c == '?' && leading > 0)
                {
                    leading = 0;
                    i++;
                }
                else
                    break;
            }
            if (leading > 0 && i < n && proForma[i] == '-')
            {
                if (leading > 1) throw Stacked(proForma, "N", leading);
                i++;
            }
            if (modified && i == n)
                throw new Tdp.ProFormaParseException("The string has a modification but no residue.");

            for (bool residue = false; i < n; i++)
            {
                char c = proForma[i];
                if (c is '[' or '{' or '<')
                {
                    int end = SkipGroup(proForma, i);
                    if (end < 0) return;
                    i = end - 1;
                    continue;
                }
                if (c != '-')
                {
                    residue |= char.IsLetter(c);
                    continue;
                }
                if (!residue)
                    throw new Tdp.ProFormaParseException($"No residue before the '-' at position {i}.");

                int j = i + 1, trailing = 0;
                while (j < n && proForma[j] == '[')
                {
                    j = SkipGroup(proForma, j);
                    if (j < 0) return;
                    trailing++;
                }
                if (trailing > 1) throw Stacked(proForma, "C", trailing);
                if (trailing == 1 && j < n && char.IsLetter(proForma[j]))
                    throw new Tdp.ProFormaParseException($"Unexpected content at position {j} after the C-terminal modification.");
                return;
            }
        }

        private static ProFormaUnsupportedException Stacked(string proForma, string side, int count) => new(
            $"'{proForma}' has {count} modifications on the {side}-terminus. That is valid ProForma 2.1 (section 6.3), but a " +
            $"ProForma term holds one modification per terminus, so reading it would lose {count - 1}.",
            $"Multiple {side}-terminal modifications ({count} found).");

        /// <summary>The index just past the group that opens at <paramref name="start"/>; -1 if it never closes.</summary>
        private static int SkipGroup(string text, int start)
        {
            int depth = 0;
            for (int k = start; k < text.Length; k++)
            {
                if (text[k] is '[' or '{' or '<') depth++;
                else if (text[k] is ']' or '}' or '>')
                {
                    depth--;
                    if (depth == 0) return k + 1;
                    if (depth < 0) return -1;
                }
            }
            return -1;
        }
    }
}
