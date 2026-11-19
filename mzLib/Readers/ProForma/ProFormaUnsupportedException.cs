using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// A ProForma string that is valid but that a <see cref="Tdp.ProFormaTerm"/> cannot hold without losing part of it.
    /// It derives from <see cref="Tdp.ProFormaParseException"/>, so a caller that catches parse failures still catches it;
    /// a caller that wants to tell "valid but unsupported" from "not valid ProForma" catches this type first.
    /// </summary>
    /// <remarks>
    /// <see cref="ProFormaReader"/> threw it for several modifications on one terminus while the SDK (0.0.299) could not
    /// hold them. Since TopDownProteomics 0.0.301 the term holds them, so the reader no longer throws it; Layer 2
    /// (<see cref="ProFormaSequenceParser"/>, <see cref="ProFormaConverter"/>) reports them as incompatible instead. The
    /// type stays so that code catching it still compiles.
    /// </remarks>
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
