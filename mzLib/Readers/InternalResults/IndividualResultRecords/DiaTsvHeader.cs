namespace Readers
{
    /// <summary>
    /// The column names MetaMorpheus's DIA tables use where no psmtsv column carries the concept. A column that does
    /// reuse a psmtsv name (<c>File Name</c>, <c>Full Sequence</c>, <c>Score</c>, ...) takes it from
    /// <see cref="SpectrumMatchFromTsvHeader"/>, so the two formats cannot drift apart. Both
    /// <see cref="DiaPrecursorFromTsv"/> and <see cref="DiaPeptideFromTsv"/> name their columns from these constants.
    /// </summary>
    public static class DiaTsvHeader
    {
        /// <summary>Precursor q-value from target-decoy competition within one run.</summary>
        public const string QValuePrecursorRun = "QValue_Precursor_Run";

        /// <summary>Precursor q-value from target-decoy competition across all runs searched together.</summary>
        public const string QValuePrecursorGlobal = "QValue_Precursor_Global";

        /// <summary>Peptide q-value from target-decoy competition among peptides, across all runs searched together.</summary>
        public const string QValuePeptideGlobal = "QValue_Peptide_Global";

        /// <summary>The library's iRT for the precursor.</summary>
        public const string LibraryIrt = "LibraryIrt";

        /// <summary>The apex's retention time on the library's iRT scale.</summary>
        public const string ApexIrt = "ApexIrt";

        /// <summary>The apex's retention time in the run, minutes.</summary>
        public const string ApexRtMin = "ApexRtMin";

        /// <summary>The one-based scan number of the apex MS2 scan.</summary>
        public const string ApexScanNumber = "ApexScanNumber";

        /// <summary>The precursor's quantity (DEF-DIA-PREC): its selected fragments' summed peak areas, MS2 only.</summary>
        public const string PrecursorQuantity = "Precursor Quantity";

        /// <summary>What the precursor quantity is (DEF-DIA-KIND): Value, BelowFloor or NotLooked.</summary>
        public const string QuantityState = "Quantity State";

        /// <summary>The floor the precursor quantity's fragments were judged against.</summary>
        public const string QuantityFloor = "Quantity Floor";

        /// <summary>The precursor's MS1 quantity (DEF-DIA-MS1).</summary>
        public const string Ms1Quantity = "MS1 Quantity";

        /// <summary>What the MS1 quantity is (DEF-DIA-KIND).</summary>
        public const string Ms1QuantityState = "MS1 Quantity State";

        /// <summary>The floor the MS1 quantity was judged against.</summary>
        public const string Ms1QuantityFloor = "MS1 Quantity Floor";
    }
}
