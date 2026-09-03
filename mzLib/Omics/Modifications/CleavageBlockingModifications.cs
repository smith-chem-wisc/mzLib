using System.Linq;
using Omics.Digestion;

namespace Omics.Modifications
{
    /// <summary>
    /// Curated classification of modifications that neutralise or mask the side-chain charge of a
    /// trypsin-family cleavage residue (Lys/Arg) enough that the protease would not have cleaved
    /// after it -- N6-succinyllysine and N6-acetyllysine being the common cases. Backs
    /// <see cref="Modification.BlocksCleavage"/>.
    ///
    /// Trypsin's specificity comes from the positively charged K/R side chain binding Asp189 in the
    /// S1 pocket; acylating the lysine epsilon-amine removes that charge (succinylation reverses it,
    /// adding a carboxylate), so the protease does not cleave. This is why acetylome and succinylome
    /// workflows are dominated by missed cleavages at the modified residue.
    ///
    /// This is a name-keyed heuristic rather than a chemistry calculation: no modification database
    /// (UniProt ptmlist, Unimod, PSI-MOD) encodes "blocks protease cleavage", so the set has to be
    /// curated somewhere. Matching is done on substrings of the modification id so that both the
    /// UniProt "N6-&lt;acyl&gt;lysine" naming and the Unimod "&lt;Acyl&gt;" short naming are covered.
    ///
    /// The methyl series (mono/di/tri-methyllysine) is deliberately EXCLUDED: methylation retains the
    /// positive charge and is sterically impaired rather than abolished, so in practice it shows up
    /// as a missed cleavage rather than an impossible one. Revisit if evidence says otherwise.
    ///
    /// Charge is the criterion, which is why glycosylation is absent from this set even though a
    /// glycan is far bulkier than an acetyl group. An N-glycan sits on N in an N-X-S/T sequon and an
    /// O-glycan on S/T -- neither is a trypsin-family cleavage residue, so there is no charge to
    /// remove and the residue gate below rejects them. Where a glycan does suppress cleavage it does
    /// so STERICALLY, usually from a neighbouring residue, and partially: the effect is a raised
    /// missed-cleavage rate, not an impossible peptidoform. That is the opposite correction to the
    /// one this class drives -- it calls for MORE missed-cleavage tolerance, not for a peptidoform to
    /// be dropped -- so modelling it here would lose real peptides. See <see cref="BlocksCleavageBy"/>
    /// for the protease half of the question.
    /// </summary>
    public static class CleavageBlockingModifications
    {
        /// <summary>
        /// Acyl (and related) groups on the lysine epsilon-amine that remove its positive charge.
        /// </summary>
        private static readonly string[] ChargeNeutralizingAcylStems =
        {
            "acetyl", "succinyl", "malonyl", "glutaryl", "crotonyl", "propionyl", "butyryl",
            "formyl", "carbamyl", "carbamoyl", "hydroxyisobutyryl", "benzoyl", "lactyl", "biotinyl",
        };

        /// <summary>
        /// The residue a modification sits on -- the one whose cleavage it would abolish -- or the null
        /// character when the modification has no target motif to read it from. A motif may carry
        /// lower-case context letters around the modified residue (an N-glycan's "Nxs"), so the
        /// modified residue is the upper-case one.
        /// </summary>
        public static char BlockedResidue(Modification modification) =>
            modification?.Target is null ? '\0' : modification.Target.Motif.FirstOrDefault(char.IsUpper);

        /// <summary>
        /// True if <paramref name="modification"/> sits on a Lys/Arg and neutralises that residue's
        /// cleavage-directing side-chain charge. False for the methyl series, and for any modification
        /// whose target is not a cleavage residue (such a modification cannot block the cleavage).
        ///
        /// This answers the CHEMISTRY question alone. Whether the configured protease actually cleaves
        /// after that residue is a separate question -- see <see cref="BlocksCleavageBy"/>, which is
        /// what digestion consults.
        /// </summary>
        public static bool NeutralizesCleavageResidue(Modification modification)
        {
            // Only Lys/Arg direct trypsin-family cleavage, so a modification anywhere else cannot block it.
            char targetResidue = BlockedResidue(modification);
            if (targetResidue != 'K' && targetResidue != 'R')
                return false;

            // Only a SIDE-CHAIN modification neutralises the cleavage-directing charge. The same acyl
            // group can appear as a protein N-terminal (alpha-amine) modification -- e.g. UniProt ships
            // both "N6-acetyllysine" (epsilon-amine, LocationRestriction "Anywhere.", genuinely blocking)
            // and "N-acetyllysine" (alpha-amine, "N-terminal.") -- and the terminal form leaves the
            // lysine side chain charged, so trypsin still cleaves. Restrict to the side-chain ("Anywhere.")
            // form so a terminal acylation is not misclassified.
            if (modification.LocationRestriction != "Anywhere.")
                return false;

            string id = (modification.OriginalId ?? modification.IdWithMotif ?? string.Empty).ToLowerInvariant();
            if (id.Length == 0)
                return false;

            // Citrullination (deimination) converts arginine's guanidinium to a neutral ureido group,
            // removing the positive charge trypsin recognises -- the one common Arg-side blocking mod.
            if (targetResidue == 'R')
                return id.Contains("citrullin") || id.Contains("deimin");

            // --- Lysine epsilon-amine chemistry below ---

            // The ubiquitin / NEDD8 / SUMO remnant left after tryptic digestion is a Gly-Gly stub,
            // which both neutralises the amine and is sterically prohibitive.
            if (id == "gg" || id.Contains("gly-gly") || id.Contains("diglycyl") || id.Contains("ubiquit"))
                return true;

            // The acyl test has the last word, and deliberately runs after nothing that could veto it.
            // An earlier revision excluded any name containing "methyl" BEFORE reaching here, which
            // misclassified the two Unimod entries carrying both chemistries on the same lysine --
            // "Methyl+Acetyl:2H(3)" and "Methyl:2H(3)+Acetyl:2H(3)", both site K, position Anywhere.
            // The acetyl on the epsilon-amine is what decides: the charge is gone whatever else the
            // residue also carries.
            //
            // The methyl series still classifies as non-blocking, and needs no guard to do it: plain
            // Methyl, Dimethyl and Trimethyl contain no acyl stem, so they fall through to false here.
            return ChargeNeutralizingAcylStems.Any(stem => id.Contains(stem));
        }

        /// <summary>
        /// True when <paramref name="modification"/> abolishes a cleavage that <paramref name="agent"/>
        /// would otherwise have performed: the chemistry blocks (see
        /// <see cref="NeutralizesCleavageResidue"/>) AND the agent actually cuts C-terminal to the
        /// residue the modification sits on. This is the predicate digestion consults, and both halves
        /// are needed -- an acetylated lysine abolishes nothing in a Glu-C digest, which never cleaved
        /// after that lysine to begin with, and counting it as blocked would discount a missed cleavage
        /// the peptide really does have.
        /// </summary>
        /// <remarks>
        /// The direction matters. Only C-terminal cleavage is recognised, because that is the only
        /// direction the digestion-side correction models: the peptidoform it drops is one whose
        /// C-TERMINUS is a blocked cut. A protease that cuts N-terminal to its recognition residue
        /// (Asp-N's "|D", Lys-N's "|K") has the mirror-image problem -- there a blocked residue
        /// invalidates the peptide's N-terminus instead -- and is left inert rather than half-corrected.
        ///
        /// Residue-level, not context-level: a protease whose motif restricts cleavage by surrounding
        /// sequence (trypsin|P's "K[P]|") still reports true for K here, so a blocked K that sits before
        /// a proline -- and was therefore never a site -- can still discount a missed cleavage. That
        /// approximation is bounded (the count is clamped at zero) and is all that remains of the much
        /// wider one this method removes.
        /// </remarks>
        public static bool BlocksCleavageBy(Modification modification, DigestionAgent agent) =>
            agent is not null
            && modification is not null
            && modification.BlocksCleavage
            && agent.CleavesCTerminalTo(BlockedResidue(modification));
    }
}
