using System.Linq;
using Omics.Digestion;

namespace Omics.Modifications
{
    /// <summary>
    /// Which glycosylation class a modification belongs to, for protease rules whose cleavage REQUIRES
    /// a glycan rather than being blocked by one. <see cref="None"/> covers every non-glycan
    /// modification, which is the overwhelming majority.
    /// </summary>
    /// <remarks>
    /// Two members rather than one flag because the distinction decides whether a rule fires at all:
    /// the O-glycoproteases (StcE, OpeRATOR/OgpA, IMPa, SmE, BT4244, AM0627, ZmpB/ZmpC, CpaA,
    /// AMUC_1438) require an O-glycan on Ser/Thr, while flavastacin requires an N-glycan on Asn. A
    /// single "is a glycan" answer would let an N-glycan satisfy an O-glycoprotease, which is wrong in
    /// both directions.
    ///
    /// <see cref="Other"/> exists because Unimod's own vocabulary has a third value ("Other
    /// glycosylation", 22 specificities in the shipped unimod.xml) covering C- and S-linked glycans,
    /// which no protease rule here keys on but which must not be silently folded into N or O.
    /// </remarks>
    public enum GlycosylationClass
    {
        None = 0,
        NLinked,
        OLinked,
        Other,
    }

    /// <summary>
    /// Curated classification of modifications that a protease may REQUIRE in order to cleave -- the
    /// mirror image of <see cref="CleavageBlockingModifications"/>. The entry point is
    /// <see cref="ClassifyGlycosylation"/>; digestion consults <see cref="Satisfies"/>.
    ///
    /// The glycoproteases invert the usual assumption. For trypsin a modification on the cleavage
    /// residue abolishes the cut; for an O-glycoprotease the cut does not exist UNLESS a glycan is
    /// present at the position the rule governs -- P1' for most of the OgpA family, P2 for StcE, both
    /// P1 and P1' for the bis-glycan enzymes. An unglycosylated Ser/Thr is simply not a site.
    ///
    /// Why a CLASS and not a modification id: the requirement is chemical, not nominal. OpeRATOR
    /// strongly prefers core 1; AMUC_1438 requires Tn only and is killed by any extended glycan; SmE
    /// accommodates sialylated core-1/core-2 and fucosylated ABO antigens; StcE takes O-GalNAc with
    /// core 1/core 2 but does NOT cleave O-GlcNAc or O-mannose. No single
    /// <see cref="Modification.IdWithMotif"/> can express any of that, so the requirement names a set.
    /// </summary>
    /// <remarks>
    /// WHAT THIS CANNOT DO YET, deliberately stated rather than approximated.
    ///
    /// The class resolved here is N-linked / O-linked / other. The finer taxonomy the rules above
    /// actually want -- core 1 vs core 2 vs Tn, O-GalNAc vs O-GlcNAc vs O-mannose -- does not exist in
    /// mzLib OR in MetaMorpheus: MetaMorpheus's GlycanType enum has exactly two members, N_glycan and
    /// O_glycan, and there is no monosaccharide-composition model on Modification at all. So a rule
    /// keyed on this class is correct but COARSE: it will admit an O-GlcNAc where StcE would not
    /// cleave. That over-admission costs search space and a little FDR burden; it does not invent
    /// identifications, and it is strictly better than today's behaviour, which ignores the glycan
    /// requirement entirely and admits every unglycosylated Ser/Thr as a site.
    ///
    /// WHERE THE CLASSIFICATION COMES FROM, in priority order. Each source is consulted because none
    /// alone is sufficient:
    ///
    /// 1. <see cref="Modification.ModificationType"/> equal to "N-linked glycosylation" or "O-linked
    ///    glycosylation". This is authoritative: it is the value MetaMorpheus's Glycan type sets on
    ///    every glycan it builds, and a glyco search's glycans are the ones that matter here.
    /// 2. <see cref="Modification.FeatureType"/> equal to "CARBOHYD" -- UniProt's feature key, on 165
    ///    shipped ptmlist entries. UniProt does not say which linkage, so the target residue decides:
    ///    Asn is N-linked, Ser/Thr O-linked.
    /// 3. Nothing else. In particular there is no id-substring fallback, unlike the blocking set's acyl
    ///    stems -- glycan names are compositional ("H2N2A2F1", "Hex(1)HexNAc(1)") rather than
    ///    descriptive, so substring matching would be guesswork with a high false-positive rate on
    ///    ordinary mods.
    ///
    /// This is thinner than it should be, and the reason is a gap in mzLib rather than in the data.
    /// The bundled unimod.xml carries a `classification` attribute on 839 specificities -- 611
    /// "O-linked glycosylation", 206 "N-linked glycosylation", 22 "Other glycosylation" -- and the
    /// generated deserializer already models it. ModificationLoader reads that file, has the value in
    /// hand, and drops it on the floor, hardcoding ModificationType to "Unimod". Carrying it through
    /// would classify every Unimod glycan here for free and make source 2 a fallback rather than a
    /// mainstay. That is a separate, additive change.
    /// </remarks>
    public static class CleavagePromotingModifications
    {
        private const string NLinkedGlycosylationModificationType = "N-linked glycosylation";
        private const string OLinkedGlycosylationModificationType = "O-linked glycosylation";
        private const string OtherGlycosylationModificationType = "Other glycosylation";

        /// <summary>
        /// UniProt's feature key for a glycosylation site. Note the spelling: UniProt writes CARBOHYD,
        /// and ptmlist ships 165 such entries.
        /// </summary>
        private const string CarbohydrateFeatureType = "CARBOHYD";

        /// <summary>
        /// The residue a modification sits on -- the one whose cleavage it would enable -- or the null
        /// character when the modification has no target motif to read it from. A motif may carry
        /// lower-case context letters around the modified residue (an N-glycan's "Nxs"), so the
        /// modified residue is the upper-case one.
        ///
        /// Identical in shape to <see cref="CleavageBlockingModifications.BlockedResidue"/>; kept
        /// separate because the two classes answer opposite questions and sharing the member would
        /// couple them for no gain.
        /// </summary>
        public static char PromotedResidue(Modification modification) =>
            modification?.Target is null ? '\0' : modification.Target.Motif.FirstOrDefault(char.IsUpper);

        /// <summary>
        /// The glycosylation class of <paramref name="modification"/>, or
        /// <see cref="GlycosylationClass.None"/> when it is not a glycan (or cannot be shown to be one
        /// from the fields available). See the remarks on <see cref="CleavagePromotingModifications"/>
        /// for the sources consulted and what is deliberately not attempted.
        ///
        /// This answers the CHEMISTRY question alone. Whether the configured protease requires a glycan
        /// at the position this one occupies is a separate question that digestion decides.
        /// </summary>
        public static GlycosylationClass ClassifyGlycosylation(Modification modification)
        {
            if (modification is null)
                return GlycosylationClass.None;

            // Source 1: the modification type, which is what MetaMorpheus's Glycan sets and is the only
            // authoritative statement of linkage available on a Modification today.
            switch (modification.ModificationType)
            {
                case NLinkedGlycosylationModificationType:
                    return GlycosylationClass.NLinked;
                case OLinkedGlycosylationModificationType:
                    return GlycosylationClass.OLinked;
                case OtherGlycosylationModificationType:
                    return GlycosylationClass.Other;
            }

            // Source 2: UniProt's CARBOHYD feature key says "this is a glycosylation site" without
            // saying which linkage, so the target residue decides. Asn carries N-linked glycans;
            // Ser/Thr carry O-linked ones. Anything else that UniProt calls CARBOHYD (C-mannosylated
            // Trp, S-linked Cys) is real but is not a linkage any protease rule here keys on, so it is
            // Other rather than guessed into N or O.
            if (modification.FeatureType == CarbohydrateFeatureType)
            {
                switch (PromotedResidue(modification))
                {
                    case 'N':
                        return GlycosylationClass.NLinked;
                    case 'S':
                    case 'T':
                        return GlycosylationClass.OLinked;
                    default:
                        return GlycosylationClass.Other;
                }
            }

            return GlycosylationClass.None;
        }

        /// <summary>
        /// True when <paramref name="modification"/> is a glycan of the class
        /// <paramref name="required"/>. <see cref="GlycosylationClass.None"/> is never satisfied by
        /// anything, so a rule that requires nothing must not be expressed as requiring None -- it must
        /// not consult this method at all.
        /// </summary>
        public static bool Satisfies(Modification modification, GlycosylationClass required) =>
            required != GlycosylationClass.None && ClassifyGlycosylation(modification) == required;

        /// <summary>
        /// True when <paramref name="modification"/> could satisfy a cleavage requirement that
        /// <paramref name="agent"/> actually has -- the chemistry matches (see <see cref="Satisfies"/>)
        /// AND at least one of the agent's motifs demands that class at some subsite. This is the
        /// predicate digestion consults, and both halves are needed: an O-glycan enables nothing in a
        /// trypsin digest, which requires no modification at all, and treating it as relevant there
        /// would buy generation slack no peptidoform can spend.
        /// </summary>
        /// <remarks>
        /// The mirror of <see cref="CleavageBlockingModifications.BlocksCleavageBy"/>, and deliberately
        /// shaped the same way: a static pair predicate taking the modification first and the agent
        /// second, matching <see cref="ModificationLocalization.ModFits"/> and the rest of
        /// <c>Omics.Modifications</c>. Null-guarded and returning false rather than throwing, because a
        /// predicate on the digestion path must not be the thing that fails a run.
        ///
        /// This answers "is this modification relevant to this agent at all", NOT "is this particular
        /// cut justified". The second question needs a position as well as a class, and is answered by
        /// <see cref="Digestion.DigestionProduct"/> when it discharges a generated peptidoform.
        /// </remarks>
        public static bool SatisfiesCleavageRequirementOf(Modification modification, DigestionAgent agent)
        {
            if (modification is null || agent?.DigestionMotifs is null)
                return false;

            foreach (DigestionMotif motif in agent.DigestionMotifs)
            {
                if (motif?.CleavageRequirement is not null && motif.CleavageRequirement.IsSatisfiedBy(modification))
                {
                    return true;
                }
            }

            return false;
        }
    }
}
