using System.Collections.Generic;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// One [Term] stanza of go.obo: the term itself and the two containment edges GO propagates over.
    ///
    /// Only is_a and part_of are kept as parents. Every other relationship -- has_part, regulates,
    /// occurs_in -- is deliberately not an edge: has_part points down the hierarchy and regulates is not
    /// containment, so following either would make an annotation imply terms it does not.
    /// </summary>
    public sealed class GeneOntologyTerm
    {
        internal GeneOntologyTerm(string id, string name, GoAspect aspect, bool isObsolete,
            IReadOnlyList<string> altIds, IReadOnlyList<string> replacedBy,
            IReadOnlyList<string> isAParents, IReadOnlyList<string> partOfParents)
        {
            Id = id;
            Name = name;
            Aspect = aspect;
            IsObsolete = isObsolete;
            AltIds = altIds;
            ReplacedBy = replacedBy;
            IsAParents = isAParents;
            PartOfParents = partOfParents;
        }

        /// <summary>The primary GO id, e.g. "GO:0005743".</summary>
        public string Id { get; }

        /// <summary>The term name, e.g. "mitochondrial inner membrane". Empty when the stanza has none.</summary>
        public string Name { get; }

        /// <summary>From the stanza's namespace; Unknown when the namespace is absent or not one of GO's three.</summary>
        public GoAspect Aspect { get; }

        /// <summary>
        /// True for an obsolete term. It stays resolvable by id, so a file citing it still reads, but it has
        /// no parents and is never an ancestor. See <see cref="ReplacedBy"/> for where GO moved it.
        /// </summary>
        public bool IsObsolete { get; }

        /// <summary>Secondary ids merged into this term. They resolve to it in <see cref="GeneOntologyGraph.TryGetTerm"/>.</summary>
        public IReadOnlyList<string> AltIds { get; }

        /// <summary>For an obsolete term, the term(s) GO says replace it. Empty otherwise.</summary>
        public IReadOnlyList<string> ReplacedBy { get; }

        /// <summary>Direct is_a parents, in file order.</summary>
        public IReadOnlyList<string> IsAParents { get; }

        /// <summary>Direct part_of parents, in file order.</summary>
        public IReadOnlyList<string> PartOfParents { get; }
    }
}
