using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// mzLib entry point for writing a <see cref="Tdp.ProFormaTerm"/> back to its canonical
    /// ProForma 2.0 string (Layer 1, term -> string). Delegates to the TopDownProteomics
    /// reference writer, which emits v2.0 canonical form.
    /// </summary>
    public static class ProFormaWriter
    {
        /// <summary>
        /// Serializes a term to its canonical ProForma 2.0 string.
        /// </summary>
        /// <param name="term">The term to serialize.</param>
        /// <returns>The canonical ProForma string.</returns>
        public static string Write(Tdp.ProFormaTerm term)
            => new Tdp.ProFormaWriter().WriteString(term);
    }
}
