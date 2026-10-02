namespace Omics.Digestion
{
    /// <summary>
    /// Where a <see cref="DigestionMotif"/>'s recognition sequence sits relative to the bond it severs.
    /// Read from <see cref="DigestionMotif.Side"/>.
    /// </summary>
    /// <remarks>
    /// THREE members, not two, and that is the whole point of the type. A motif's sidedness is encoded
    /// in <see cref="DigestionMotif.CutIndex"/> -- the position of the '|' within the de-bracketed motif
    /// string -- and that has three regimes, not the two a naive reading assumes:
    ///
    /// <list type="bullet">
    /// <item>CutIndex == 0 -- the motif lies wholly after the cut. Asp-N "|D", Lys-N "|K",
    /// CNBr_N "|M", RNase_MC1 "|U". <see cref="NTerminal"/>.</item>
    /// <item>CutIndex == InducingCleavage.Length -- the motif lies wholly before the cut. Trypsin "K|",
    /// Lys-C "K|", Glu-C "E|". <see cref="CTerminal"/>.</item>
    /// <item>0 &lt; CutIndex &lt; InducingCleavage.Length -- the motif STRADDLES the bond, naming
    /// residues on both sides. Collagenase "GPX|GPX" (CutIndex 3), StcE-trypsin "TX|T" (CutIndex 2),
    /// colicin_E5 "G|U" (CutIndex 1). <see cref="Straddling"/>.</item>
    /// </list>
    ///
    /// The third case is not hypothetical and not rare enough to approximate: all four motifs named
    /// above ship in proteases.tsv and rnases.tsv today. Code that classifies sidedness by testing
    /// CutIndex against a literal -- "CutIndex == 1 means C-terminal", "CutIndex == 0 means
    /// N-terminal" -- silently misses every straddling motif, which is why this distinction is a type
    /// rather than a comparison written out at each call site.
    ///
    /// Note that a straddling motif is still a C-terminal CUTTER in the sense
    /// <see cref="DigestionMotif.CleavesCTerminalTo"/> asks about: it severs the bond after some
    /// residue it names. Side describes where the recognition sequence lies; CleavesCTerminalTo asks
    /// which residue the cut falls after. Both are true of "TX|T", and neither implies the other.
    /// </remarks>
    public enum CleavageSide
    {
        /// <summary>
        /// The cut falls before every residue the motif names, so the motif describes only the prime
        /// (P1', P2', ...) side. An N-terminal cutter: Asp-N "|D".
        /// </summary>
        NTerminal,

        /// <summary>
        /// The cut falls after every residue the motif names, so the motif describes only the
        /// non-prime (P1, P2, ...) side. A C-terminal cutter: trypsin "K|".
        /// </summary>
        CTerminal,

        /// <summary>
        /// The cut falls inside the motif, which therefore names residues on both sides of the bond.
        /// Collagenase "GPX|GPX", StcE-trypsin "TX|T", colicin_E5 "G|U".
        /// </summary>
        Straddling,
    }
}
