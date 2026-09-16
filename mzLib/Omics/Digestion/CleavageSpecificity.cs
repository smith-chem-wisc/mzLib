namespace Omics.Digestion
{
    /// <summary>
    /// How many termini of a digestion product are made by the digestion agent's cleavage rules.
    /// </summary>
    /// <remarks>
    /// <para>The same enum is used in three places, and a value does not mean quite the same thing in each:</para>
    /// <list type="number">
    /// <item><description><b>A digestion agent's own specificity</b> (<see cref="DigestionAgent.CleavageSpecificity"/>, the
    /// Specificity column of proteases.tsv): Full, Semi, None (no cleavage at all: top-down, peptidomics), SingleN or
    /// SingleC.</description></item>
    /// <item><description><b>The kind of search a caller asks for</b> (<see cref="IDigestionParams.SearchModeType"/>):
    /// Full, Semi or None. Together with <see cref="IDigestionParams.FragmentationTerminus"/> it decides whether digestion
    /// returns peptides or seeds for an engine that trims them after the search; see the remarks on
    /// <c>DigestionParams.SearchModeType</c>.</description></item>
    /// <item><description><b>The label on one digestion product</b>
    /// (<see cref="DigestionProduct.CleavageSpecificityForFdrCategory"/>): Full or Semi for peptides, SingleN or SingleC
    /// for non-specific seeds, and Unknown for a peptide a search engine trimmed from a seed, whose specificity is worked
    /// out afterwards. MetaMorpheus uses this label to split PSMs into separate FDR categories.</description></item>
    /// </list>
    /// </remarks>
    public enum CleavageSpecificity
    {
        /// <summary>
        /// No terminus needs to be made by the cleavage rules. As an agent's specificity: no cleavage (top-down). As a
        /// search mode: a non-specific search, which in digestion means singleN or singleC seeds, never a list of
        /// non-specific peptides.
        /// </summary>
        None,

        /// <summary>
        /// At least one terminus made by the cleavage rules (the other may be anywhere); includes the fully specific
        /// products. As a search mode it means semi-specific peptides with FragmentationTerminus Both, and seeds with N or C.
        /// </summary>
        Semi,

        /// <summary>Both termini made by the cleavage rules. The default search mode.</summary>
        Full,

        /// <summary>
        /// The singleN agent, and its products: seeds with a fixed N-terminus at every position, whose C-terminal end a
        /// non-specific search engine decides afterwards.
        /// </summary>
        SingleN,

        /// <summary>
        /// The singleC agent, and its products: seeds with a fixed C-terminus at every position, whose N-terminal end a
        /// non-specific search engine decides afterwards.
        /// </summary>
        SingleC,

        /// <summary>
        /// The label of a peptide a search engine trimmed from a seed after the search (fast semi-specific and non-specific
        /// searches); its real specificity is determined from its termini afterwards.
        /// </summary>
        Unknown
    }
}
