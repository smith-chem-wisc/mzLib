using Omics.Modifications;

namespace Omics.Digestion
{
    /// <summary>
    /// A modification a <see cref="DigestionMotif"/> REQUIRES at one Schechter-Berger subsite before it
    /// will sever its bond -- the mirror image of <see cref="Modifications.CleavageBlockingModifications"/>,
    /// which describes a modification that abolishes a cleavage that would otherwise happen.
    ///
    /// The O-glycoproteases invert the usual assumption. For trypsin a modification on the cleavage
    /// residue removes a site that the bare sequence already had; for StcE or OpeRATOR the site does not
    /// exist at all unless a glycan is present, so an unglycosylated Ser/Thr is simply not a site. A
    /// sequence motif cannot say that, which is why <c>proteases.tsv</c> records that its StcE entries
    /// "OVER-DIGEST relative to the real enzyme".
    /// </summary>
    /// <remarks>
    /// <para><b>Why a subsite address and not a flag.</b> The required glycan does not sit in the same
    /// place for every enzyme, and the differences are not cosmetic:</para>
    /// <list type="bullet">
    /// <item><description>OpeRATOR/OgpA, IMPa, SmE, BT4244: <c>X | T/S*</c> -- the glycan is at
    /// <b>P1'</b>, the residue immediately after the bond.</description></item>
    /// <item><description>StcE: <c>T/S* X | T/S</c> -- the glycan is at <b>P2</b>, two residues BEFORE
    /// the bond, and the motif straddles it. Note the P1' residue carries no asterisk: a glycan there
    /// is permitted and never required.</description></item>
    /// <item><description>AM0627: <c>T/S* | T/S*</c> -- glycans on BOTH sides.</description></item>
    /// <item><description>ZmpC: <c>T/S* X X X | X</c> -- the glycan is at <b>P4</b>, four residues
    /// upstream. Nothing attached to the residues flanking the bond can express that.</description></item>
    /// </list>
    /// <para>So the address is (side, position), where position is an ordinary Schechter-Berger subsite
    /// number and is NOT limited to 1. Use <see cref="NonPrime"/> and <see cref="Prime"/> rather than
    /// constructing one directly, so the call site reads as the enzymology does:
    /// <c>CleavageRequirement.NonPrime(2, GlycosylationClass.OLinked)</c> is StcE.</para>
    ///
    /// <para><b>Do not confuse the side here with <see cref="CleavageSide"/>.</b> That enum says where a
    /// motif's recognition sequence sits relative to the bond and has THREE values, because a motif may
    /// straddle. A subsite is a single residue and is genuinely either non-prime or prime, so two values
    /// are exhaustive here and no third case is being overlooked.</para>
    ///
    /// <para><b>What this deliberately cannot express.</b> The requirement names a glycosylation CLASS,
    /// not a structure. Real enzymes discriminate far more finely -- OpeRATOR needs core 1 and is refused
    /// by the Tn antigen, AMUC_1438 needs Tn and is refused by anything larger, ZmpB requires an
    /// alpha-2,6-sialyl branch -- but neither mzLib nor MetaMorpheus models monosaccharide composition on
    /// a <see cref="Modification"/>, so a structure-level rule cannot be written yet. A class-level rule
    /// is coarse in one direction only: it admits glycoforms the real enzyme would refuse. That costs
    /// search space, never an identification, and it is strictly better than today, where the glycan is
    /// ignored entirely and every bare Ser/Thr is treated as a site.</para>
    /// </remarks>
    public sealed class CleavageRequirement
    {
        private CleavageRequirement(bool isPrimeSide, int subsite, GlycosylationClass requiredClass)
        {
            IsPrimeSide = isPrimeSide;
            Subsite = subsite;
            RequiredClass = requiredClass;
        }

        /// <summary>
        /// True when the constrained residue lies AFTER the severed bond (P1', P2', ...), false when it
        /// lies before it (P1, P2, ...). Which side it is decides which of the two digestion products
        /// can answer for the cut: see <see cref="DigestionProduct"/>'s discharge.
        /// </summary>
        public bool IsPrimeSide { get; }

        /// <summary>
        /// The Schechter-Berger subsite number, counting outward from the bond starting at 1. P1 and P1'
        /// are the residues either side of it. Not capped at 1 -- ZmpC needs 4.
        /// </summary>
        public int Subsite { get; }

        /// <summary>The glycosylation class the residue at that subsite must carry.</summary>
        public GlycosylationClass RequiredClass { get; }

        /// <summary>
        /// A requirement on the NON-PRIME side, before the bond. <c>NonPrime(2, OLinked)</c> is StcE's:
        /// an O-glycan two residues before the cut.
        /// </summary>
        public static CleavageRequirement NonPrime(int subsite, GlycosylationClass requiredClass) =>
            new(isPrimeSide: false, subsite, requiredClass);

        /// <summary>
        /// A requirement on the PRIME side, after the bond. <c>Prime(1, OLinked)</c> is the OgpA family's:
        /// an O-glycan on the residue the cut exposes as a new N-terminus.
        /// </summary>
        public static CleavageRequirement Prime(int subsite, GlycosylationClass requiredClass) =>
            new(isPrimeSide: true, subsite, requiredClass);

        /// <summary>
        /// True when <paramref name="modification"/> is a glycan of the class this requirement demands.
        /// A modification whose <see cref="Modification.ModificationType"/> is "Protease" can never
        /// satisfy it, whatever its class: those are the cleavage modifications a protease leaves on the
        /// terminus it just created, and <see cref="DigestionProduct"/> places them at exactly the keys a
        /// P1 or P1' check reads -- key 2 for the N-terminal one and key length+1 for the C-terminal one.
        /// Without this guard a cleavage modification would be read as the required glycan and justify
        /// its own cut.
        /// </summary>
        public bool IsSatisfiedBy(Modification modification) =>
            modification is not null
            && modification.ModificationType != ProteaseModificationType
            && CleavagePromotingModifications.Satisfies(modification, RequiredClass);

        /// <summary>
        /// The <see cref="Modification.ModificationType"/> that marks a protease-associated cleavage
        /// modification. Matches the literal <see cref="DigestionProduct"/> tests when it places them.
        /// </summary>
        internal const string ProteaseModificationType = "Protease";

        public override string ToString() =>
            "P" + Subsite + (IsPrimeSide ? "'" : string.Empty) + " requires " + RequiredClass;
    }
}
