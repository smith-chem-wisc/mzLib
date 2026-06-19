using System.Collections.Generic;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// A single peptidoform that lives above the single ProForma term: one or more chains joined by
    /// <c>//</c> (inter-chain crosslinks and branches), with an optional <c>/charge</c> and ion adducts.
    /// Each chain is a <see cref="Tdp.ProFormaTerm"/> handled by the wrapped single-term SDK parser/writer.
    /// </summary>
    public class ProFormaPeptidoform
    {
        /// <summary>The chains of this peptidoform, in order; joined by <c>//</c> when more than one.</summary>
        public IReadOnlyList<Tdp.ProFormaTerm> Chains { get; }

        /// <summary>The charge state (signed), or <c>null</c> when the peptidoform carries no <c>/z</c>.</summary>
        public int? Charge { get; }

        /// <summary>
        /// Raw ion-adduct content between the brackets following the charge (for example
        /// <c>+2Na+,+H+</c> in <c>/2[+2Na+,+H+]</c>), or <c>null</c> when none is present.
        /// </summary>
        public string IonAdducts { get; }

        public ProFormaPeptidoform(IReadOnlyList<Tdp.ProFormaTerm> chains, int? charge, string ionAdducts)
        {
            Chains = chains;
            Charge = charge;
            IonAdducts = ionAdducts;
        }
    }

    /// <summary>
    /// A full ProForma 2.0 string at the level above a single proteoform term: one or more
    /// <em>chimeric</em> peptidoforms joined by <c>+</c> (HUPO-PSI ProForma 2.0 §7.2), each a charge-bearing
    /// (§7.1) set of crosslinked / branched chains joined by <c>//</c> (§4.2.3.2, §4.2.4). This is the
    /// "layer above the single-term model"; per-chain parsing and writing are delegated to the wrapped
    /// <c>TopDownProteomics</c> single-term implementation.
    /// </summary>
    public class ProFormaMultiTerm
    {
        /// <summary>The chimeric peptidoforms, in order; joined by <c>+</c> when more than one.</summary>
        public IReadOnlyList<ProFormaPeptidoform> Peptidoforms { get; }

        public ProFormaMultiTerm(IReadOnlyList<ProFormaPeptidoform> peptidoforms)
        {
            Peptidoforms = peptidoforms;
        }
    }
}
