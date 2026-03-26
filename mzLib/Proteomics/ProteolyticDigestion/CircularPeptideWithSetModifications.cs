using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;

namespace Proteomics.ProteolyticDigestion
{
    /// <summary>
    /// Represents a circular (head-to-tail cyclized) peptide with set modifications,
    /// always bound to a <see cref="CircularProtein"/> parent.
    ///
    /// NUMBERING
    /// ---------
    /// All residue numbering is inherited from the parent <see cref="CircularProtein"/>.
    /// Residue 1 is the first character of <see cref="CircularProtein.CircularSequence"/>
    /// (the canonical origin — the lexicographically smallest rotation of the ring).
    /// Fragment annotations [s, e] use this numbering, with e > N for wrap-around fragments.
    ///
    /// MASS
    /// ----
    /// A circular peptide has no free termini. Its monoisotopic mass equals the sum of
    /// residue masses plus modification masses, with no added H₂O:
    ///   MonoisotopicMass = base.MonoisotopicMass − H₂O
    ///
    /// FRAGMENTATION
    /// -------------
    /// Only <see cref="FragmentInternally"/> is supported. A single backbone cleavage
    /// opens the ring and produces a linear ion mass-equivalent to the precursor — no
    /// additional sequence information is gained and these ions are not scored.
    ///
    /// <see cref="FragmentInternally"/> enumerates all contiguous sub-sequences of length
    /// [minLength, peptideLength−1] within this peptide's span in the parent ring,
    /// including wrap-around fragments for full-ring peptides. Residue masses are looked
    /// up from the parent ring via a doubled prefix-sum for O(1) sub-sequence mass
    /// calculation.
    /// </summary>
    [Serializable]
    public class CircularPeptideWithSetModifications : PeptideWithSetModifications
    {
        // ── Mass constant ─────────────────────────────────────────────────────
        private static readonly double WaterMonoisotopicMass =
            PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2
            + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        // ── Cached mass arrays (lazy-initialized, immutable after first use) ──
        private double[] _cachedRingMasses;
        private double[] _cachedDoubledPrefix;

        // ── Parent protein (typed) ────────────────────────────────────────────

        /// <summary>
        /// The parent <see cref="CircularProtein"/>. Always non-null; the constructor
        /// enforces this. Use this reference to access the canonical numbering system
        /// and ring-level properties (e.g., <see cref="CircularProtein.CyclicMonoisotopicMass"/>).
        /// </summary>
        public CircularProtein CircularParent { get; }

        // ── Constructor ───────────────────────────────────────────────────────

        /// <summary>
        /// Constructs a <see cref="CircularPeptideWithSetModifications"/> from a
        /// <see cref="CircularProtein"/> parent.
        ///
        /// No canonicalization is performed here: the parent <see cref="CircularProtein"/>
        /// guarantees its <see cref="Protein.BaseSequence"/> and all sub-sequences derived
        /// from <see cref="CircularProtein.Digest"/> are already in canonical form.
        ///
        /// The intended way to obtain instances is via <see cref="CircularProtein.Digest"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">
        /// Thrown if <paramref name="protein"/> is null.
        /// </exception>
        /// <exception cref="ArgumentException">
        /// Thrown if <paramref name="protein"/> is not a <see cref="CircularProtein"/>.
        /// </exception>
        public CircularPeptideWithSetModifications(
            Protein protein,
            IDigestionParams digestionParams,
            int oneBasedStartResidueInProtein,
            int oneBasedEndResidueInProtein,
            CleavageSpecificity cleavageSpecificity,
            string peptideDescription,
            int missedCleavages,
            Dictionary<int, Modification> allModsOneIsNterminus,
            int numFixedMods,
            string baseSequence = null)
            : base(ValidateProtein(protein), digestionParams,
                   oneBasedStartResidueInProtein, oneBasedEndResidueInProtein,
                   cleavageSpecificity, peptideDescription, missedCleavages,
                   allModsOneIsNterminus, numFixedMods, baseSequence)
        {
            // ValidateProtein already guarantees protein is a non-null CircularProtein.
            CircularParent = (CircularProtein)protein;

            // Terminal modifications are chemically invalid on a circular peptide
            // (no free N- or C-terminus). Reject them early so that mass accounting
            // in MonoisotopicMass and FragmentInternally stays consistent.
            if (allModsOneIsNterminus != null)
            {
                int peptideLen = oneBasedEndResidueInProtein - oneBasedStartResidueInProtein + 1;
                if (allModsOneIsNterminus.ContainsKey(1))
                    throw new ArgumentException(
                        "N-terminal modifications are not valid for circular peptides (no free N-terminus).",
                        nameof(allModsOneIsNterminus));
                if (allModsOneIsNterminus.ContainsKey(peptideLen + 2))
                    throw new ArgumentException(
                        "C-terminal modifications are not valid for circular peptides (no free C-terminus).",
                        nameof(allModsOneIsNterminus));
            }
        }

        /// <summary>
        /// Validates that <paramref name="protein"/> is a non-null <see cref="CircularProtein"/>
        /// before the base constructor executes. Called inline in the constructor initializer
        /// so that callers receive descriptive exceptions instead of opaque crashes from the
        /// base class.
        /// </summary>
        private static Protein ValidateProtein(Protein protein)
        {
            if (protein is null)
                throw new ArgumentNullException(nameof(protein),
                    $"{nameof(CircularPeptideWithSetModifications)} requires a non-null protein.");

            if (protein is not CircularProtein)
                throw new ArgumentException(
                    $"{nameof(CircularPeptideWithSetModifications)} requires a " +
                    $"{nameof(CircularProtein)} parent, but received {protein.GetType().Name}.",
                    nameof(protein));

            return protein;
        }

        // ── Mass override ─────────────────────────────────────────────────────

        /// <summary>
        /// Monoisotopic mass of the circular peptide: sum of residue and modification
        /// masses with no added water (no free termini).
        /// Equals the parent protein's <see cref="CircularProtein.CyclicMonoisotopicMass"/>
        /// when unmodified.
        /// </summary>
        /// <summary>
        /// Monoisotopic mass of the circular peptide: sum of residue and modification
        /// masses with no added water (no free termini).
        /// Overrides the base implementation to subtract H2O, reflecting the absence of
        /// free termini in a head-to-tail cyclized peptide.
        /// </summary>
        public override double MonoisotopicMass =>
            (double)ClassExtensions.RoundedDouble(
                base.MonoisotopicMass - WaterMonoisotopicMass);

        // ── Fragmentation ─────────────────────────────────────────────────────

        /// <summary>
        /// Generates internal fragment ions for this circular peptide.
        ///
        /// NUMBERING
        /// ---------
        /// Fragment positions are expressed in the parent <see cref="CircularProtein"/>'s
        /// canonical 1-based numbering. For a fragment of length L starting at canonical
        /// residue s:
        ///   FragmentNumber           = s          (1-based canonical start)
        ///   SecondaryFragmentNumber  = s + L − 1  (end; exceeds N for wrap-around)
        ///   ResiduePosition          = L          (fragment length)
        ///
        /// MASS
        /// ----
        ///   neutralMass = sum(residueMasses in parent ring over [s, s+L-1])
        ///                 + nTermCap + cTermCap − H₂O
        ///
        /// SPAN RESTRICTION
        /// ----------------
        /// Only fragments within this peptide's span in the ring are generated — both
        /// cleavage sites must fall within the peptide. For a full-ring peptide
        /// (Length == N) this covers all ring positions including wrap-around. For a
        /// sub-peptide (Length &lt; N) wrap-around internal fragments are excluded because
        /// they would require a cleavage outside the peptide's span.
        /// </summary>
        public override void FragmentInternally(
            DissociationType dissociationType,
            int minLengthOfFragments,
            List<Product> products,
            IFragmentationParams? fragmentationParams = null)
        {
            if (minLengthOfFragments < 1)
                throw new ArgumentOutOfRangeException(nameof(minLengthOfFragments),
                    minLengthOfFragments,
                    $"{nameof(minLengthOfFragments)} must be at least 1.");

            products.Clear();

            int peptideLength = BaseSequence.Length;
            int ringLength = CircularParent.Length;

            // Internal fragments require at least 2 residues: one cleavage on each side
            // of the fragment within the peptide span.
            if (minLengthOfFragments < 2 || minLengthOfFragments >= peptideLength)
                return;

            // ── Build doubled prefix sum over the parent ring ─────────────────
            // Masses come from the parent ring (not just the peptide subsequence)
            // so that mod positions and residue masses are correct for the full ring.
            // Arrays are cached because this peptide's state is immutable after
            // construction, so the values never change between calls.
            double[] ringMasses = GetOrBuildRingMasses();
            int peptideStartInRing = OneBasedStartResidueInProtein - 1; // 0-based

            double[] doubledPrefix = GetOrBuildDoubledPrefix(ringMasses, peptideStartInRing, peptideLength, ringLength);

            // ── Ion type caps ─────────────────────────────────────────────────
            var massCaps = DissociationTypeCollection
                .GetNAndCTerminalMassShiftsForDissociationType(dissociationType);

            List<ProductType> nTermProductTypes =
                DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(
                    dissociationType, FragmentationTerminus.N);

            List<ProductType> cTermProductTypes =
                DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(
                    dissociationType, FragmentationTerminus.C);

            // ── Enumerate internal fragments ──────────────────────────────────
            // localStart: 0-based offset from this peptide's first residue.
            //
            // Sub-peptide (peptideLength < ringLength):
            //   Both cleavage sites must fall within [0, peptideLength), so
            //   maxLength = peptideLength - localStart - 1, preventing wrap-around.
            //
            // Full-ring peptide (peptideLength == ringLength):
            //   Every pair of backbone cleavage sites on the ring is valid. A fragment
            //   may wrap around the origin, so maxLength = peptideLength - 1 for all
            //   starting positions (the doubled prefix-sum handles the wrap).

            bool isFullRing = peptideLength == ringLength;

            for (int localStart = 0; localStart < peptideLength; localStart++)
            {
                int oneBasedStart = OneBasedStartResidueInProtein + localStart;
                int maxLength = isFullRing
                    ? peptideLength - 1
                    : peptideLength - localStart - 1;

                for (int length = minLengthOfFragments; length <= maxLength; length++)
                {
                    double fragmentResidueMass =
                        doubledPrefix[localStart + length] - doubledPrefix[localStart];

                    int oneBasedEnd = oneBasedStart + length - 1;

                    for (int i = 0; i < nTermProductTypes.Count; i++)
                    {
                        double nTermCap = massCaps.Item1[i];
                        for (int j = 0; j < cTermProductTypes.Count; j++)
                        {
                            double cTermCap = massCaps.Item2[j];

                            products.Add(new Product(
                                productType: cTermProductTypes[j],
                                terminus: FragmentationTerminus.None,
                                neutralMass: fragmentResidueMass + nTermCap + cTermCap - WaterMonoisotopicMass,
                                fragmentNumber: oneBasedStart,
                                residuePosition: length,
                                neutralLoss: 0,
                                secondaryProductType: nTermProductTypes[i],
                                secondaryFragmentNumber: oneBasedEnd));
                        }
                    }
                }
            }
        }

        // ── Private helper ────────────────────────────────────────────────────

        /// <summary>
        /// Builds a residue mass array of length N (the parent ring length) from
        /// the parent ring's canonical sequence, with this peptide's modifications
        /// mapped back to their correct ring positions.
        ///
        /// This is the single mass array used by <see cref="FragmentInternally"/>.
        /// Using the full ring (not just the peptide subsequence) ensures that:
        ///   (a) residue masses are always read from the correct ring position, and
        ///   (b) the doubled prefix-sum correctly handles wrap-around for full-ring
        ///       peptides without re-entering the peptide's own residues out of order.
        ///
        /// Modifications in <see cref="PeptideWithSetModifications.AllModsOneIsNterminus"/>
        /// use keys relative to the peptide's local sequence. We convert each to a
        /// 0-based ring index:
        ///   ringIndex = (OneBasedStartResidueInProtein − 1 + localIndex) % ringLength
        /// </summary>
        /// <summary>
        /// Returns the cached ring mass array, building it on first access.
        /// </summary>
        private double[] GetOrBuildRingMasses()
        {
            return _cachedRingMasses ??= BuildParentRingMassArray();
        }

        /// <summary>
        /// Returns the cached doubled prefix-sum array, building it on first access.
        /// </summary>
        private double[] GetOrBuildDoubledPrefix(double[] ringMasses, int peptideStartInRing, int peptideLength, int ringLength)
        {
            if (_cachedDoubledPrefix != null)
                return _cachedDoubledPrefix;

            var doubledPrefix = new double[2 * peptideLength + 1];
            for (int i = 0; i < 2 * peptideLength; i++)
            {
                int ringIndex = (peptideStartInRing + i) % ringLength;
                doubledPrefix[i + 1] = doubledPrefix[i] + ringMasses[ringIndex];
            }

            _cachedDoubledPrefix = doubledPrefix;
            return _cachedDoubledPrefix;
        }

        private double[] BuildParentRingMassArray()
        {
            int ringLength = CircularParent.Length;
            string ringSeq = CircularParent.BaseSequence;
            double[] masses = new double[ringLength];

            for (int i = 0; i < ringLength; i++)
            {
                if (!Residue.TryGetResidue(ringSeq[i], out Residue res))
                    throw new InvalidOperationException(
                        $"Unrecognized amino acid '{ringSeq[i]}' at position {i + 1} " +
                        $"in the ring sequence of protein '{CircularParent.Accession}'. " +
                        $"Cannot compute fragment masses.");

                masses[i] = res.MonoisotopicMass;
            }

            // Map each mod from its peptide-local key to its ring position.
            // Keys:  1       → local index 0  (N-term mod)
            //        i+2     → local index i  (side-chain, 0-based)
            //        n+2     → local index n-1 (C-term mod, n = peptide length)
            int n = BaseSequence.Length;
            int startInRing = OneBasedStartResidueInProtein - 1;

            foreach (var kvp in AllModsOneIsNterminus)
            {
                // Only side-chain mods (keys 2..n+1) are valid; terminal mods
                // are rejected by the constructor.
                int localIndex = kvp.Key - 2;

                if (localIndex < 0 || localIndex >= n)
                    throw new InvalidOperationException(
                        $"Unexpected modification key {kvp.Key} for circular peptide of length {n}. " +
                        "Only side-chain modifications (keys 2 to n+1) are supported.");

                int ringIndex = (startInRing + localIndex) % ringLength;
                masses[ringIndex] += kvp.Value.MonoisotopicMass.Value;
            }

            return masses;
        }

        // ── ToString ──────────────────────────────────────────────────────────

        public override string ToString() =>
            $"[Circular] {FullSequence} [{OneBasedStartResidueInProtein}-{OneBasedEndResidueInProtein}]";
    }
}