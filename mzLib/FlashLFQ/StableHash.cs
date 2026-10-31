using System;
using System.Buffers;
using System.Text;

namespace FlashLFQ
{
    /// <summary>
    /// A hash that gives the same value for the same text in every run, process and machine. It is FNV-1a
    /// (64-bit) over the text's UTF-8 bytes, followed by MurmurHash3's 64-bit finalizer so that texts that
    /// differ only near their end still spread over the whole range.
    ///
    /// <see cref="string.GetHashCode()"/> cannot be used where a result must be reproducible: .NET seeds it
    /// randomly in every process.
    /// </summary>
    internal static class StableHash
    {
        internal const ulong FnvOffsetBasis = 14695981039346656037UL;
        private const ulong FnvPrime = 1099511628211UL;
        private const byte FieldSeparator = 0x1F;

        /// <summary>
        /// Continues an FNV-1a state over the UTF-8 bytes of <paramref name="text"/>.
        /// </summary>
        internal static ulong Append(ulong state, string text)
        {
            if (string.IsNullOrEmpty(text))
            {
                return state;
            }

            int maxBytes = Encoding.UTF8.GetMaxByteCount(text.Length);
            byte[] rented = null;
            Span<byte> buffer = maxBytes <= 512 ? stackalloc byte[maxBytes] : (rented = ArrayPool<byte>.Shared.Rent(maxBytes));
            try
            {
                int length = Encoding.UTF8.GetBytes(text, buffer);
                for (int i = 0; i < length; i++)
                {
                    state = (state ^ buffer[i]) * FnvPrime;
                }
                return state;
            }
            finally
            {
                if (rented != null)
                {
                    ArrayPool<byte>.Shared.Return(rented);
                }
            }
        }

        /// <summary>
        /// Continues an FNV-1a state over one separator byte, so that ("AB", "C") and ("A", "BC") hash apart.
        /// </summary>
        internal static ulong AppendSeparator(ulong state) => (state ^ FieldSeparator) * FnvPrime;

        /// <summary>
        /// FNV-1a (64-bit) of the UTF-8 bytes of <paramref name="text"/>, without the finalizer.
        /// </summary>
        internal static ulong Fnv1a64(string text) => Append(FnvOffsetBasis, text);

        /// <summary>
        /// MurmurHash3's 64-bit finalizer (fmix64).
        /// </summary>
        internal static ulong Mix(ulong k)
        {
            k ^= k >> 33;
            k *= 0xff51afd7ed558ccdUL;
            k ^= k >> 33;
            k *= 0xc4ceb9fe1a85ec53UL;
            k ^= k >> 33;
            return k;
        }
    }
}
