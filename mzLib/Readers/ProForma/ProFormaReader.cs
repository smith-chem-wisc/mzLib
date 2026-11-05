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
    /// Since TopDownProteomics 0.0.301 the term holds several modifications on one terminus (ProForma 2.1, section 6.3)
    /// in <see cref="Tdp.ProFormaTerm.NTerminalModifications"/> and <see cref="Tdp.ProFormaTerm.CTerminalModifications"/>,
    /// and every labile group in <see cref="Tdp.ProFormaTerm.LabileModifications"/>. The SDK also refuses a residue after
    /// the C-terminal modification and a string with no residue. 0.0.299 lost or moved those modifications without an
    /// error, which is why this reader used to refuse those shapes itself. <see cref="Tdp.ProFormaTerm.NTerminalDescriptors"/>
    /// and <see cref="Tdp.ProFormaTerm.CTerminalDescriptors"/> are only the first modification of each terminus: a caller
    /// that can hold one per terminus must check the counts (Layer 2 does).
    /// </remarks>
    public static class ProFormaReader
    {
        /// <summary>
        /// Parses a ProForma string into its term representation.
        /// </summary>
        /// <param name="proFormaString">A ProForma 2.0 string.</param>
        /// <returns>The parsed <see cref="Tdp.ProFormaTerm"/>.</returns>
        /// <exception cref="Tdp.ProFormaParseException">Thrown when the input is not valid ProForma, including a residue
        /// after the C-terminal modification and a string with no residue.</exception>
        public static Tdp.ProFormaTerm Read(string proFormaString)
            => new Tdp.ProFormaParser().ParseString(proFormaString);

        /// <summary>
        /// The first terminus of <paramref name="term"/> that carries more than one modification (ProForma 2.1, section 6.3),
        /// as "N" or "C" with its count; null when each terminus has at most one. For callers that hold one modification
        /// per terminus, which must refuse such a term rather than read only <see cref="Tdp.ProFormaTerm.NTerminalDescriptors"/>
        /// or <see cref="Tdp.ProFormaTerm.CTerminalDescriptors"/>.
        /// </summary>
        internal static (string Side, int Count)? StackedTerminus(Tdp.ProFormaTerm term)
        {
            int n = term.NTerminalModifications?.Count ?? 0;
            if (n > 1) return ("N", n);
            int c = term.CTerminalModifications?.Count ?? 0;
            return c > 1 ? ("C", c) : null;
        }
    }
}
