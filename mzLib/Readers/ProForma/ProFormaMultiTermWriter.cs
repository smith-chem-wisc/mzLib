using System.Globalization;
using System.Linq;
using System.Text;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// mzLib entry point for writing a <see cref="ProFormaMultiTerm"/> back to its canonical ProForma 2.0
    /// string. Chains are serialized by the wrapped single-term <see cref="ProFormaWriter"/> and rejoined
    /// with the inter-chain <c>//</c>, charge <c>/z[adducts]</c>, and chimeric <c>+</c> operators.
    /// The single-term canonical form is the SDK's, so this writer is canonically idempotent wherever the
    /// SDK writer is.
    /// </summary>
    public static class ProFormaMultiTermWriter
    {
        public static string Write(ProFormaMultiTerm multiTerm)
        {
            return string.Join("+", multiTerm.Peptidoforms.Select(WritePeptidoform));
        }

        private static string WritePeptidoform(ProFormaPeptidoform peptidoform)
        {
            var sb = new StringBuilder();
            sb.Append(string.Join("//", peptidoform.Chains.Select(ProFormaWriter.Write)));

            if (peptidoform.Charge.HasValue)
            {
                sb.Append('/').Append(peptidoform.Charge.Value.ToString(CultureInfo.InvariantCulture));
                if (peptidoform.IonAdducts != null)
                    sb.Append('[').Append(peptidoform.IonAdducts).Append(']');
            }
            return sb.ToString();
        }
    }
}
