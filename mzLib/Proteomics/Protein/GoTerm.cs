using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    /// <summary>
    /// The three Gene Ontology aspects. UniProt writes these as a single-character prefix on the
    /// dbReference's "term" property -- "C:cytoplasm", "F:NAD binding", "P:glycolytic process".
    /// </summary>
    public enum GoAspect
    {
        /// <summary>
        /// No aspect prefix was present, or it was not one UniProt uses. Deliberately distinct from
        /// the three real aspects rather than guessed at: an organelle claim built on an inferred
        /// aspect is worse than one that admits it does not know.
        /// </summary>
        Unknown = 0,
        BiologicalProcess,
        CellularComponent,
        MolecularFunction
    }

    /// <summary>
    /// One Gene Ontology annotation on a protein, deduplicated by GO id.
    ///
    /// This is a derived view over Protein.DatabaseReferences, not a stored field -- GO already
    /// arrives there because ProteinXmlEntry keeps every dbReference generically. See
    /// Protein.GoTerms.
    /// </summary>
    public class GoTerm
    {
        public GoTerm(string id, GoAspect aspect, string name,
            IEnumerable<string> evidenceCodes = null, IEnumerable<string> projects = null)
        {
            Id = id ?? "";
            Aspect = aspect;
            Name = name ?? "";
            EvidenceCodes = new HashSet<string>(evidenceCodes ?? Enumerable.Empty<string>(), StringComparer.Ordinal);
            Projects = new HashSet<string>(projects ?? Enumerable.Empty<string>(), StringComparer.Ordinal);
        }

        /// <summary>
        /// The GO accession, e.g. "GO:0005737". Unique within a protein's GoTerms.
        /// </summary>
        public string Id { get; }

        /// <summary>
        /// Cellular component, molecular function or biological process -- or Unknown when the
        /// database did not say.
        /// </summary>
        public GoAspect Aspect { get; }

        /// <summary>
        /// The term name with its aspect prefix removed, e.g. "cytoplasm" for "C:cytoplasm".
        /// Empty when the reference carried no term property.
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Every ECO evidence code supporting this term on this protein. A SET, because UniProt
        /// repeats one GO id once per line of evidence and the codes are what a downstream filter
        /// (drop everything supported only by IEA-like codes) acts on. Empty is normal.
        /// </summary>
        public IReadOnlyCollection<string> EvidenceCodes { get; }

        /// <summary>
        /// The annotating projects (UniProtKB, MGI, HPA, ...), unioned across the same repeats.
        /// </summary>
        public IReadOnlyCollection<string> Projects { get; }

        /// <summary>
        /// Projects a flat list of database references onto GO terms: keeps only Type == "GO",
        /// deduplicates by id, and unions the evidence codes and projects of the repeats.
        ///
        /// Filtering by type is mandatory, not defensive -- the same list also holds PubMed, DOI,
        /// EMBL and organism references, and the XML parser applies no depth guard, so references
        /// nested inside other elements arrive here too.
        /// </summary>
        public static IReadOnlyList<GoTerm> FromDatabaseReferences(IEnumerable<DatabaseReference> references)
        {
            if (references == null)
            {
                return Array.Empty<GoTerm>();
            }

            return references
                .Where(r => r != null && r.Type == Protein.GeneOntologyDatabaseReferenceType)
                .GroupBy(r => r.Id ?? "", StringComparer.Ordinal)
                .Select(BuildTerm)
                .ToList();
        }

        private static GoTerm BuildTerm(IGrouping<string, DatabaseReference> sameId)
        {
            var evidence = new List<string>();
            var projects = new List<string>();
            string rawTerm = null;

            foreach (var property in sameId.SelectMany(PropertiesOf))
            {
                // Matched by type, never by position: the properties are an unordered bag, and
                // ProteinDbWriter re-sorts them on write, so a round-tripped file already presents
                // them in a different order from the one UniProt shipped.
                switch (property.Item1)
                {
                    case "evidence":
                        AddIfPresent(evidence, property.Item2);
                        break;
                    case "project":
                        AddIfPresent(projects, property.Item2);
                        break;
                    case "term":
                        rawTerm = rawTerm ?? property.Item2;
                        break;
                }
            }

            SplitTerm(rawTerm, out GoAspect aspect, out string name);
            return new GoTerm(sameId.Key, aspect, name, evidence, projects);
        }

        private static IEnumerable<Tuple<string, string>> PropertiesOf(DatabaseReference reference) =>
            reference.Properties ?? Enumerable.Empty<Tuple<string, string>>();

        private static void AddIfPresent(List<string> into, string value)
        {
            if (!string.IsNullOrEmpty(value))
            {
                into.Add(value);
            }
        }

        /// <summary>
        /// Splits "C:cytoplasm" into its aspect and its name. Only a single-character C/F/P prefix
        /// counts: term names legitimately contain colons, so anything else is left whole and the
        /// aspect is reported Unknown rather than invented.
        /// </summary>
        private static void SplitTerm(string rawTerm, out GoAspect aspect, out string name)
        {
            aspect = GoAspect.Unknown;
            name = rawTerm ?? "";

            if (name.Length < 2 || name[1] != ':')
            {
                return;
            }

            switch (name[0])
            {
                case 'C': aspect = GoAspect.CellularComponent; break;
                case 'F': aspect = GoAspect.MolecularFunction; break;
                case 'P': aspect = GoAspect.BiologicalProcess; break;
                default: return; // unrecognised prefix: keep the whole string as the name
            }

            name = name.Substring(2);
        }
    }
}
