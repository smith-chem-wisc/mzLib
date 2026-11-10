using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    /// <summary>
    /// An O-glycan a motif needs, or must not have, at one subsite before it cuts, written as the literature does: StcE is
    /// "P2:O-glycan", and IMPa, which does not cut between two adjacent glycosites, adds "P1:!O-glycan".
    /// </summary>
    public sealed class CleavageRequirement
    {
        public CleavageRequirement(int subsite, bool isPrimeSide, bool isForbidden = false)
        {
            Subsite = subsite;
            IsPrimeSide = isPrimeSide;
            IsForbidden = isForbidden;
        }

        /// <summary>Counted outward from the cut: P1 is the residue before it, P1' the residue after it.</summary>
        public int Subsite { get; }

        public bool IsPrimeSide { get; }

        /// <summary>True when an O-glycan at the subsite prevents the cut instead of being needed for it.</summary>
        public bool IsForbidden { get; }

        public int ConstrainedResidue(int cutAfterOneBasedResidue) =>
            IsPrimeSide ? cutAfterOneBasedResidue + Subsite : cutAfterOneBasedResidue - Subsite + 1;

        public bool IsSatisfiedBy(Modification modification) =>
            modification.ModificationType == "O-linked glycosylation"
            || (modification.FeatureType == "CARBOHYD" && modification.Target?.Motif.FirstOrDefault(char.IsUpper) is 'S' or 'T');

        /// <summary>Reads a "Cleavage Requirement" cell of proteases.tsv: subsites separated by ";", e.g. "P1':O-glycan;P1:!O-glycan".</summary>
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
                    || !parts[1].TrimStart('!').Equals("O-glycan", StringComparison.OrdinalIgnoreCase))
                {
                    throw new MzLibException($"Unrecognized cleavage requirement '{clause}'. Expected a subsite and O-glycan, for example P2:O-glycan, or P1:!O-glycan where a glycan prevents the cut.");
                }
                requirements.Add(new CleavageRequirement(number, parts[0].EndsWith('\''), parts[1].StartsWith('!')));
            }
            return requirements;
        }
    }
}
