using System.Linq;

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
        /// True if <paramref name="modification"/> sits on a Lys/Arg and neutralises that residue's
        /// cleavage-directing side-chain charge. False for the methyl series, and for any modification
        /// whose target is not a cleavage residue (such a modification cannot block the cleavage).
        /// </summary>
        public static bool NeutralizesCleavageResidue(Modification modification)
        {
            if (modification?.Target is null)
                return false;

            // Only Lys/Arg direct trypsin-family cleavage, so a modification anywhere else cannot block it.
            char targetResidue = modification.Target.Motif.FirstOrDefault(char.IsUpper);
            if (targetResidue != 'K' && targetResidue != 'R')
                return false;

            string id = (modification.OriginalId ?? modification.IdWithMotif ?? string.Empty).ToLowerInvariant();
            if (id.Length == 0)
                return false;

            // Methylation retains the charge -- excluded by design (see class remarks).
            if (id.Contains("methyl"))
                return false;

            // The ubiquitin / NEDD8 / SUMO remnant left after tryptic digestion is a Gly-Gly stub,
            // which both neutralises the amine and is sterically prohibitive.
            if (id == "gg" || id.Contains("gly-gly") || id.Contains("diglycyl") || id.Contains("ubiquit"))
                return true;

            return ChargeNeutralizingAcylStems.Any(stem => id.Contains(stem));
        }
    }
}
