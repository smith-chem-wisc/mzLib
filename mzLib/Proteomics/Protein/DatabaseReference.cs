using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class DatabaseReference
    {
        /// <summary>
        /// DatabaseReference constructor, which takes the type and ID strings of the reference, and a list of properties. Each property contains the "type" and "value" of the property as Item1 and Item2 of the Tuple.
        /// </summary>
        /// <param name="type"></param>
        /// <param name="id"></param>
        /// <param name="properties"></param>
        public DatabaseReference(string type, string id, IEnumerable<Tuple<string, string>> properties)
        {
            Type = type ?? "";
            Id = id ?? "";
            Properties = properties ?? new List<Tuple<string, string>>();
        }

        /// <summary>
        /// A reference that UniProt scopes to one isoform, as in
        /// &lt;dbReference type="Reactome" id="R-HSA-5250924"&gt;&lt;molecule id="O00159-3"/&gt;...&lt;/dbReference&gt;.
        /// </summary>
        /// <param name="type"></param>
        /// <param name="id"></param>
        /// <param name="properties"></param>
        /// <param name="molecule">The isoform accession the reference is scoped to, or null (or empty) when it covers
        /// the whole entry.</param>
        public DatabaseReference(string type, string id, IEnumerable<Tuple<string, string>> properties, string molecule)
            : this(type, id, properties)
        {
            Molecule = string.IsNullOrEmpty(molecule) ? null : molecule;
        }

        /// <summary>
        /// The isoform this reference is scoped to (UniProt's &lt;molecule id&gt;, e.g. "O00159-3"), or null when it
        /// covers the whole entry. An isoform-scoped reference does not apply to the entry's other isoforms, so it must
        /// never be read as the entry's: in the UniProt 2026_03 human proteome, 610 Reactome pathways are listed for
        /// one isoform only.
        /// </summary>
        public string Molecule { get; }

        /// <summary>
        /// dbRef type, e.g. "GO" for GO terms
        /// </summary>
        public string Type { get; }

        /// <summary>
        /// dbRef ID string
        /// </summary>
        public string Id { get; }

        /// <summary>
        /// Each database reference contains a list of properties. Item1 of this Tuple is the "type", and Item2 is the "value" of the property.
        /// </summary>
        public IEnumerable<Tuple<string, string>> Properties { get; }

        public override bool Equals(object obj)
        {
            DatabaseReference d = obj as DatabaseReference;
            return d != null
                && (d.Type == null && Type == null || d.Type.Equals(Type))
                && (d.Id == null && Id == null || d.Id.Equals(Id))
                && string.Equals(d.Molecule, Molecule, StringComparison.Ordinal)
                && d.Properties.OrderBy(x => x).SequenceEqual(Properties.OrderBy(x => x));
        }

        public override int GetHashCode()
        {
            return Type.GetHashCode() ^ Id.GetHashCode(); // null handled in constructor
        }
    }
}