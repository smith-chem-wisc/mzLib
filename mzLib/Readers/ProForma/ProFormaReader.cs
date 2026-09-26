using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// mzLib entry point for parsing a HUPO-PSI ProForma 2.0 proteoform string into a
    /// <see cref="Tdp.ProFormaTerm"/> (Layer 1, lossless string -> term).
    /// Parsing is delegated to the TopDownProteomics reference implementation; this facade
    /// exists so mzLib/MetaMorpheus callers depend on a stable mzLib type, not the SDK directly.
    /// </summary>
    public static class ProFormaReader
    {
        /// <summary>
        /// Parses a ProForma string into its term representation.
        /// </summary>
        /// <param name="proFormaString">A ProForma 2.0 string.</param>
        /// <returns>The parsed <see cref="Tdp.ProFormaTerm"/>.</returns>
        /// <exception cref="Tdp.ProFormaParseException">Thrown when the input is not valid ProForma.</exception>
        public static Tdp.ProFormaTerm Read(string proFormaString)
            => new Tdp.ProFormaParser().ParseString(proFormaString);
    }
}
