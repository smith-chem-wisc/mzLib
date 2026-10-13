using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// A ProForma string that is valid but that a <see cref="Tdp.ProFormaTerm"/> cannot hold without losing part of it,
    /// such as several modifications on one terminus (ProForma 2.1, section 6.3). It derives from
    /// <see cref="Tdp.ProFormaParseException"/>, so a caller that catches parse failures still catches it; a caller that
    /// wants to tell "valid but unsupported" from "not valid ProForma" catches this type first.
    /// </summary>
    public sealed class ProFormaUnsupportedException : Tdp.ProFormaParseException
    {
        public ProFormaUnsupportedException(string message, string incompatibleItem) : base(message)
        {
            IncompatibleItem = incompatibleItem;
        }

        /// <summary>A short statement of what could not be held, e.g. "Multiple C-terminal modifications (2 found)."</summary>
        public string IncompatibleItem { get; }
    }
}
