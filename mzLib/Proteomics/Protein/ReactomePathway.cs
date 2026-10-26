using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    /// <summary>
    /// One Reactome pathway a protein entry is a member of, as UniProt lists it.
    ///
    /// This is a derived view over Protein.DatabaseReferences, not a stored field -- UniProt's
    /// &lt;dbReference type="Reactome" id="R-HSA-70171"&gt; entries already arrive there because ProteinXmlEntry
    /// keeps every dbReference generically, each with its &lt;molecule&gt; when UniProt scopes it to one isoform.
    /// See Protein.ReactomePathways.
    ///
    /// UniProt lists the lowest-level pathways a protein takes part in, not their ancestors, so this is direct
    /// membership only: rolling up to broader pathways needs Reactome's own hierarchy file.
    /// </summary>
    public class ReactomePathway
    {
        /// <summary>The property type UniProt uses for the pathway's name on a Reactome dbReference.</summary>
        public const string PathwayNamePropertyType = "pathway name";

        public ReactomePathway(string id, string name, IEnumerable<string> molecules = null)
        {
            Id = id ?? "";
            Name = name ?? "";
            Molecules = new HashSet<string>(molecules ?? Enumerable.Empty<string>(), StringComparer.Ordinal);
        }

        /// <summary>The Reactome stable id, e.g. "R-HSA-70171". Unique within a protein's ReactomePathways.</summary>
        public string Id { get; }

        /// <summary>
        /// The pathway name as UniProt wrote it, e.g. "Glycolysis" -- UniProt's snapshot of Reactome's name.
        /// Empty when the reference carried none.
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// The isoforms the membership is scoped to (UniProt's &lt;molecule id&gt;, e.g. "P04406-2"). EMPTY means the
        /// whole entry. Not empty means the pathway is documented for those isoforms only, and does not apply to the
        /// entry's other isoforms.
        /// </summary>
        public IReadOnlyCollection<string> Molecules { get; }

        /// <summary>
        /// Projects a flat list of database references onto Reactome pathways: keeps only Type == "Reactome", one pathway
        /// per id, and unions the isoforms of repeats.
        ///
        /// Filtering by type is mandatory, not defensive -- the same list holds GO, PubMed, organism and every other
        /// reference, and the XML parser applies no depth guard. A reference without a name is still a membership and is
        /// kept with an empty name; dropping it would read as "not annotated". A repeat that covers the whole entry
        /// covers every isoform, so the merged pathway is whole-entry.
        /// </summary>
        public static IReadOnlyList<ReactomePathway> FromDatabaseReferences(IEnumerable<DatabaseReference> references)
        {
            if (references == null)
            {
                return Array.Empty<ReactomePathway>();
            }

            return references
                .Where(r => r != null && r.Type == Protein.ReactomeDatabaseReferenceType)
                .GroupBy(r => r.Id ?? "", StringComparer.Ordinal)
                .Select(BuildPathway)
                .ToList();
        }

        private static ReactomePathway BuildPathway(IGrouping<string, DatabaseReference> sameId)
        {
            // Matched by type, never by position: ProteinDbWriter re-sorts properties on write.
            string name = sameId.Select(r => PropertyValue(r, PathwayNamePropertyType)).FirstOrDefault(n => !string.IsNullOrEmpty(n));
            var molecules = sameId.Any(r => r.Molecule == null) ? Enumerable.Empty<string>() : sameId.Select(r => r.Molecule);
            return new ReactomePathway(sameId.Key, name, molecules);
        }

        private static string PropertyValue(DatabaseReference reference, string type) =>
            (reference.Properties ?? Enumerable.Empty<Tuple<string, string>>())
                .FirstOrDefault(p => p != null && p.Item1 == type)?.Item2;
    }
}
