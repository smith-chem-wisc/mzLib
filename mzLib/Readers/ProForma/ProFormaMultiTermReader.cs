using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// mzLib entry point for parsing a ProForma 2.0 string that lives <em>above</em> a single proteoform
    /// term — chimeric peptidoforms (<c>+</c>), charge states (<c>/z[adducts]</c>), and multi-chain /
    /// branch constructs (<c>//</c>). The string is split on these top-level operators (respecting
    /// bracket nesting) and each resulting chain is parsed by the wrapped single-term
    /// <see cref="ProFormaReader"/>. A plain single-term string parses to one peptidoform of one chain.
    /// </summary>
    public static class ProFormaMultiTermReader
    {
        /// <summary>Parses a (possibly multi-term) ProForma string into its <see cref="ProFormaMultiTerm"/>.</summary>
        /// <exception cref="Tdp.ProFormaParseException">If a chain is not valid ProForma.</exception>
        public static ProFormaMultiTerm Read(string proFormaString)
        {
            if (proFormaString == null) throw new ArgumentNullException(nameof(proFormaString));

            var peptidoforms = new List<ProFormaPeptidoform>();
            foreach (var pf in SplitTopLevel(proFormaString, '+'))
                peptidoforms.Add(ParsePeptidoform(pf));
            return new ProFormaMultiTerm(peptidoforms);
        }

        private static ProFormaPeptidoform ParsePeptidoform(string pf)
        {
            int chargeSlash = FindChargeSlash(pf);
            string chainsPart = chargeSlash < 0 ? pf : pf.Substring(0, chargeSlash);

            int? charge = null;
            string adducts = null;
            if (chargeSlash >= 0)
                (charge, adducts) = ParseCharge(pf.Substring(chargeSlash + 1));

            var chains = new List<Tdp.ProFormaTerm>();
            foreach (var chain in SplitDoubleSlash(chainsPart))
                chains.Add(ProFormaReader.Read(chain));

            return new ProFormaPeptidoform(chains, charge, adducts);
        }

        /// <summary>Parses the text after the charge <c>/</c>: a signed integer plus an optional <c>[adducts]</c>.</summary>
        private static (int charge, string adducts) ParseCharge(string s)
        {
            int bracket = s.IndexOf('[');
            string number = bracket < 0 ? s : s.Substring(0, bracket);
            string adducts = null;
            if (bracket >= 0)
            {
                int close = s.LastIndexOf(']');
                if (close <= bracket)
                    throw new Tdp.ProFormaParseException("Unterminated charge adduct bracket in '" + s + "'.");
                adducts = s.Substring(bracket + 1, close - bracket - 1);
            }
            if (!int.TryParse(number, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int charge))
                throw new Tdp.ProFormaParseException("Invalid charge value '" + number + "'.");
            return (charge, adducts);
        }

        // ---- bracket-depth-aware splitting --------------------------------------------------------

        private static bool IsOpen(char c) => c == '[' || c == '(' || c == '{' || c == '<';
        private static bool IsClose(char c) => c == ']' || c == ')' || c == '}' || c == '>';

        /// <summary>Splits on a single delimiter that occurs at bracket depth 0.</summary>
        private static IEnumerable<string> SplitTopLevel(string s, char delim)
        {
            var parts = new List<string>();
            var sb = new StringBuilder();
            int depth = 0;
            foreach (char c in s)
            {
                if (IsOpen(c)) depth++;
                else if (IsClose(c)) depth--;

                if (depth == 0 && c == delim)
                {
                    parts.Add(sb.ToString());
                    sb.Clear();
                }
                else sb.Append(c);
            }
            parts.Add(sb.ToString());
            return parts;
        }

        /// <summary>Splits on <c>//</c> occurring at bracket depth 0 (the inter-chain / branch separator).</summary>
        private static IEnumerable<string> SplitDoubleSlash(string s)
        {
            var parts = new List<string>();
            var sb = new StringBuilder();
            int depth = 0;
            for (int i = 0; i < s.Length; i++)
            {
                char c = s[i];
                if (IsOpen(c)) depth++;
                else if (IsClose(c)) depth--;

                if (depth == 0 && c == '/' && i + 1 < s.Length && s[i + 1] == '/')
                {
                    parts.Add(sb.ToString());
                    sb.Clear();
                    i++; // consume the second '/'
                }
                else sb.Append(c);
            }
            parts.Add(sb.ToString());
            return parts;
        }

        /// <summary>
        /// Returns the index of the depth-0 charge slash — a single <c>/</c> not part of a <c>//</c> —
        /// or -1 if the peptidoform carries no charge. Chain-separator <c>//</c> pairs are skipped.
        /// </summary>
        private static int FindChargeSlash(string s)
        {
            int depth = 0;
            for (int i = 0; i < s.Length; i++)
            {
                char c = s[i];
                if (IsOpen(c)) depth++;
                else if (IsClose(c)) depth--;
                else if (depth == 0 && c == '/')
                {
                    if (i + 1 < s.Length && s[i + 1] == '/') { i++; continue; } // chain separator
                    return i;
                }
            }
            return -1;
        }
    }
}
