#nullable enable
using System;
using System.Collections.Generic;

namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>
    /// The names of taxonomic ranks and their order from the root down. A rank is a plain string, kept exactly as the
    /// source wrote it on each <see cref="TaxonomyNode"/>; this class is how two of them are compared.
    /// </summary>
    public static class TaxonomicRank
    {
        /// <summary>The top level of cellular life (Bacteria, Archaea, Eukaryota). NCBI formerly called it "superkingdom".</summary>
        public const string Domain = "domain";
        public const string Kingdom = "kingdom";
        public const string Phylum = "phylum";
        public const string Class = "class";
        public const string Order = "order";
        public const string Family = "family";
        public const string Genus = "genus";
        public const string Species = "species";

        /// <summary>
        /// The eight ranks most analyses group by, root first: domain, kingdom, phylum, class, order, family, genus,
        /// species.
        /// </summary>
        public static readonly IReadOnlyList<string> MajorRanks = new[]
        {
            Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species
        };

        // Depth from the root; a larger number is a narrower group. Gaps leave room for the ranks between.
        private static readonly Dictionary<string, int> Depth = new(StringComparer.OrdinalIgnoreCase)
        {
            ["acellular root"] = -10,
            ["cellular root"] = -10,
            ["realm"] = -5,
            ["subrealm"] = -4,
            ["domain"] = 0,
            ["superkingdom"] = 0,
            ["kingdom"] = 10,
            ["subkingdom"] = 15,
            ["superphylum"] = 18,
            ["phylum"] = 20,
            ["subphylum"] = 25,
            ["infraphylum"] = 27,
            ["superclass"] = 28,
            ["class"] = 30,
            ["subclass"] = 35,
            ["infraclass"] = 37,
            ["cohort"] = 38,
            ["superorder"] = 39,
            ["order"] = 40,
            ["suborder"] = 45,
            ["infraorder"] = 46,
            ["parvorder"] = 47,
            ["superfamily"] = 48,
            ["family"] = 50,
            ["subfamily"] = 55,
            ["tribe"] = 58,
            ["subtribe"] = 59,
            ["genus"] = 60,
            ["subgenus"] = 65,
            ["section"] = 66,
            ["series"] = 67,
            ["species group"] = 68,
            ["species subgroup"] = 69,
            ["species"] = 70,
            ["subspecies"] = 75,
            ["varietas"] = 76,
            ["forma"] = 77,
            ["strain"] = 80,
            ["serotype"] = 80,
            ["serogroup"] = 80,
            ["biotype"] = 80,
            ["genotype"] = 80,
            ["isolate"] = 80,
        };

        /// <summary>
        /// Whether <paramref name="rank"/> is a named level with a place in the hierarchy. False for null, empty,
        /// "no rank", "clade" and any name this class does not know: those say nothing about how deep a taxon sits.
        /// </summary>
        public static bool IsRanked(string? rank) => rank != null && Depth.ContainsKey(rank.Trim());

        /// <summary>
        /// Whether two rank names mean the same level. Case and surrounding space are ignored, and "superkingdom" is
        /// the same level as "domain". Two unranked names ("clade", "no rank") are never the same level.
        /// </summary>
        public static bool AreSameRank(string? a, string? b) =>
            a != null && b != null
            && Depth.TryGetValue(a.Trim(), out int depthOfA)
            && Depth.TryGetValue(b.Trim(), out int depthOfB)
            && depthOfA == depthOfB
            && (depthOfA != 80 || string.Equals(a.Trim(), b.Trim(), StringComparison.OrdinalIgnoreCase));

        /// <summary>
        /// Orders two rank names from the root down: negative when <paramref name="a"/> is the broader group, positive
        /// when it is the narrower, zero when they are the same level. Null when either is unranked, because an
        /// unranked taxon can sit at any depth.
        /// </summary>
        public static int? Compare(string? a, string? b)
        {
            if (a == null || b == null
                || !Depth.TryGetValue(a.Trim(), out int depthOfA)
                || !Depth.TryGetValue(b.Trim(), out int depthOfB))
            {
                return null;
            }
            return depthOfA.CompareTo(depthOfB);
        }
    }
}
