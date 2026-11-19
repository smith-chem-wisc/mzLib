using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    /// <summary>
    /// An O-glycan a motif needs at one subsite before it cuts, written as the literature does: StcE is "P2:O-glycan".
    /// </summary>
    public sealed class CleavageRequirement
    {
        public CleavageRequirement(int subsite, bool isPrimeSide)
        {
            Subsite = subsite;
            IsPrimeSide = isPrimeSide;
        }

        /// <summary>Counted outward from the cut: P1 is the residue before it, P1' the residue after it.</summary>
        public int Subsite { get; }

        public bool IsPrimeSide { get; }

        public int ConstrainedResidue(int cutAfterOneBasedResidue) =>
            IsPrimeSide ? cutAfterOneBasedResidue + Subsite : cutAfterOneBasedResidue - Subsite + 1;

        public bool IsSatisfiedBy(Modification modification) =>
            modification.ModificationType == "O-linked glycosylation"
            || (modification.FeatureType == "CARBOHYD" && modification.Target?.Motif.FirstOrDefault(char.IsUpper) is 'S' or 'T');

        /// <summary>Reads a "Cleavage Requirement" cell of proteases.tsv: subsites separated by ";", e.g. "P2:O-glycan".</summary>
        public static List<CleavageRequirement> ParseAll(string text)
        {
            var requirements = new List<CleavageRequirement>();
            foreach (string clause in text.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
            {
                string[] parts = clause.Replace(" ", string.Empty).Split(':');
                string subsite = parts[0].TrimEnd('\'');
                if (parts.Length != 2
                    || !subsite.StartsWith("P", StringComparison.OrdinalIgnoreCase)
                    || !int.TryParse(subsite[1..], out int number) || number < 1
                    || !parts[1].Equals("O-glycan", StringComparison.OrdinalIgnoreCase))
                {
                    throw new MzLibException($"Unrecognized cleavage requirement '{clause}'. Expected a subsite and O-glycan, for example P2:O-glycan or P1':O-glycan.");
                }
                requirements.Add(new CleavageRequirement(number, parts[0].EndsWith('\'')));
            }
            return requirements;
        }
    }
}
