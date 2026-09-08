namespace Readers
{
    /// <summary>The kind of inconsistency found across a set of documents.</summary>
    public enum SdrfDriftKind
    {
        /// <summary>
        /// One accession, written with different names. The accession is the identity and the name
        /// is only a label, so this does not break a join -- but it means the pooled table displays
        /// the same thing two ways, and it is the earliest visible sign that two annotations were
        /// produced by different hands or different vocabulary versions.
        /// </summary>
        AccessionNameConflict,

        /// <summary>
        /// One name, mapped to different accessions. The serious one: the same word has been given
        /// two machine identities, so anything keyed on the accession silently splits, and anything
        /// keyed on the name silently merges two different things.
        /// </summary>
        NameAccessionConflict,

        /// <summary>
        /// A column annotated with controlled-vocabulary terms in some documents and with free text
        /// in others. The two do not compare equal, so half the corpus drops out of any query that
        /// assumes one or the other.
        /// </summary>
        MixedTermAndFreeText,

        /// <summary>
        /// Column names that differ only by case or spacing. SDRF comparison is case- and
        /// space-sensitive, so these are different columns as far as any consumer is concerned,
        /// and the values in them will never be compared to one another.
        /// </summary>
        ColumnNameVariant,

        /// <summary>
        /// Free-text values in one column that differ only by case or surrounding whitespace, e.g.
        /// "label free sample" and "Label Free Sample". Equality-based grouping treats them as two
        /// populations.
        /// </summary>
        ValueCaseVariant
    }

    /// <summary>One way a concept was written, and how widely.</summary>
    public sealed record SdrfDriftVariant(string Value, int Occurrences, IReadOnlyList<string> Documents)
    {
        public override string ToString() =>
            $"'{Value}' x{Occurrences} in {Documents.Count} document(s)";
    }

    /// <summary>
    /// One concept written more than one way across a collection.
    /// </summary>
    /// <param name="Kind">See <see cref="SdrfDriftKind"/>.</param>
    /// <param name="Concept">
    /// What was written inconsistently: an accession, a term name, or a column name, depending on
    /// the kind.
    /// </param>
    /// <param name="Column">The column involved, or null when the finding is about a column name itself.</param>
    /// <param name="Variants">Every spelling seen, most frequent first. The first is the corpus's own majority usage.</param>
    public sealed record SdrfDriftFinding(
        SdrfDriftKind Kind,
        string Concept,
        string? Column,
        IReadOnlyList<SdrfDriftVariant> Variants)
    {
        /// <summary>
        /// The most widely used spelling in the collection.
        ///
        /// This describes what the documents DO; it is not advice about what to write. Pooled with
        /// public data it will frequently name the wrong thing -- in the curated corpus,
        /// characteristics[organism] is free text in 963 documents against a controlled-vocabulary
        /// term in 275, so following the majority would mean abandoning the accession. Authored
        /// values come from the pinned vocabulary regardless of what the majority does.
        ///
        /// It is still the right signal for a genuine spelling drift -- casing of a term name,
        /// prefix style on an accession -- where every variant means the same thing and only one
        /// is conventional.
        /// </summary>
        public SdrfDriftVariant Majority => Variants[0];

        public override string ToString() =>
            $"{Kind} '{Concept}'" + (Column is null ? "" : $" in [{Column}]") +
            ": " + string.Join(" vs ", Variants.Select(v => v.ToString()));
    }

    /// <summary>
    /// Finds concepts that a set of SDRF documents annotated inconsistently.
    ///
    /// This is the check the whole project turns on. Every document can be individually valid --
    /// <see cref="SdrfValidator"/> can pass all of them -- and the pooled table still be unusable,
    /// because validity is a property of one file and comparability is a property of the
    /// relationship between files. Nothing in the specification, and nothing in any single-file
    /// validator, looks across documents.
    ///
    /// It has a second use. Run over the community's curated corpus it reports how everyone else
    /// writes a given concept -- which is a description of the corpus, NOT a recommendation. Where
    /// variants are merely spellings of one thing (term-name casing, accession prefix style), the
    /// majority is worth matching. Where they differ in kind, it is not: most public documents
    /// write organism as free text, and copying that would throw away the accession. Authored
    /// values are resolved from the pinned vocabulary irrespective of the majority; tolerating
    /// other people's free text is this layer's job, not the writer's.
    /// </summary>
    public static class SdrfDriftLint
    {
        /// <summary>
        /// Columns whose values are legitimately unique per row, so variation in them is data rather
        /// than drift. Comparing them would bury the real findings in noise.
        /// </summary>
        private static readonly HashSet<string> PerRowColumns = new(StringComparer.Ordinal)
        {
            "source name",
            "assay name",
            "comment[data file]",
            "comment[file uri]",
            SdrfCollection.SourceDocumentColumn
        };

        public static IReadOnlyList<SdrfDriftFinding> Analyze(SdrfCollection collection)
        {
            if (collection is null) throw new ArgumentNullException(nameof(collection));

            var findings = new List<SdrfDriftFinding>();
            findings.AddRange(FindColumnNameVariants(collection));
            findings.AddRange(FindValueDrift(collection));

            // Sorted, because the findings are built by walking Dictionary instances and .NET
            // guarantees no enumeration order. Returning them unordered is harmless for a console
            // report and fatal the moment the report becomes an artefact you diff between runs to
            // see whether the corpus is drifting -- which is the standing check this type exists to
            // support. Order by impact first so the useful entries stay at the top.
            return findings
                .OrderByDescending(f => f.Variants.Skip(1).Sum(v => v.Occurrences))
                .ThenBy(f => f.Kind)
                .ThenBy(f => f.Column ?? "", StringComparer.Ordinal)
                .ThenBy(f => f.Concept, StringComparer.Ordinal)
                .ToList();
        }

        /// <summary>
        /// Column names that collapse onto one another once case and internal spacing are ignored.
        /// </summary>
        private static IEnumerable<SdrfDriftFinding> FindColumnNameVariants(SdrfCollection collection)
        {
            // normalised name -> exact spelling -> documents that used it
            var byNormalized = new Dictionary<string, Dictionary<string, List<string>>>(StringComparer.Ordinal);

            for (int d = 0; d < collection.Count; d++)
            {
                foreach (var name in collection[d].Header.Distinct(StringComparer.Ordinal))
                {
                    if (string.IsNullOrWhiteSpace(name)) continue;
                    string key = Normalize(name);
                    if (!byNormalized.TryGetValue(key, out var spellings))
                        byNormalized[key] = spellings = new Dictionary<string, List<string>>(StringComparer.Ordinal);
                    if (!spellings.TryGetValue(name, out var documents))
                        spellings[name] = documents = new List<string>();
                    documents.Add(collection.Labels[d]);
                }
            }

            foreach (var kv in byNormalized.Where(k => k.Value.Count > 1))
            {
                var variants = kv.Value
                    .Select(s => new SdrfDriftVariant(s.Key, s.Value.Count, s.Value))
                    .OrderByDescending(v => v.Occurrences)
                    .ThenBy(v => v.Value, StringComparer.Ordinal)
                    .ToList();

                yield return new SdrfDriftFinding(
                    SdrfDriftKind.ColumnNameVariant, kv.Key, null, variants);
            }
        }

        private static IEnumerable<SdrfDriftFinding> FindValueDrift(SdrfCollection collection)
        {
            // column -> accession -> name -> documents
            var namesByAccession = new Dictionary<string, Dictionary<string, Dictionary<string, List<string>>>>(StringComparer.Ordinal);
            // column -> name -> accession -> documents
            var accessionsByName = new Dictionary<string, Dictionary<string, Dictionary<string, List<string>>>>(StringComparer.Ordinal);
            // column -> normalised free text -> exact text -> documents
            var freeTextVariants = new Dictionary<string, Dictionary<string, Dictionary<string, List<string>>>>(StringComparer.Ordinal);
            // column -> documents using terms / documents using free text
            var termDocuments = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
            var freeTextDocuments = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);

            for (int d = 0; d < collection.Count; d++)
            {
                var document = collection[d];
                string label = collection.Labels[d];
                var header = document.Header;

                foreach (var row in document.Results)
                {
                    for (int c = 0; c < row.Cells.Count && c < header.Count; c++)
                    {
                        string column = header[c];
                        string cell = row.Cells[c];

                        if (string.IsNullOrWhiteSpace(column) || string.IsNullOrWhiteSpace(cell)) continue;
                        if (PerRowColumns.Contains(column)) continue;
                        if (IsReserved(cell)) continue;

                        if (SdrfCell.TryParseTerm(cell, out var term))
                        {
                            Track(termDocuments, column, label);

                            if (!string.IsNullOrEmpty(term.Accession) && !string.IsNullOrEmpty(term.Name))
                            {
                                Add(namesByAccession, column, term.Accession, term.Name, label);
                                Add(accessionsByName, column, term.Name, term.Accession, label);
                            }
                        }
                        else
                        {
                            Track(freeTextDocuments, column, label);
                            Add(freeTextVariants, column, Normalize(cell), cell, label);
                        }
                    }
                }
            }

            foreach (var finding in Conflicts(namesByAccession, SdrfDriftKind.AccessionNameConflict))
                yield return finding;

            foreach (var finding in Conflicts(accessionsByName, SdrfDriftKind.NameAccessionConflict))
                yield return finding;

            foreach (var finding in Conflicts(freeTextVariants, SdrfDriftKind.ValueCaseVariant))
                yield return finding;

            foreach (var column in termDocuments.Keys)
            {
                if (!freeTextDocuments.TryGetValue(column, out var freeText)) continue;
                var terms = termDocuments[column];
                if (terms.Count == 0 || freeText.Count == 0) continue;

                yield return new SdrfDriftFinding(
                    SdrfDriftKind.MixedTermAndFreeText, column, column,
                    new[]
                    {
                        new SdrfDriftVariant("controlled vocabulary term", terms.Count, terms.OrderBy(x => x, StringComparer.Ordinal).ToList()),
                        new SdrfDriftVariant("free text", freeText.Count, freeText.OrderBy(x => x, StringComparer.Ordinal).ToList())
                    }
                    .OrderByDescending(v => v.Occurrences)
                    .ToList());
            }
        }

        private static IEnumerable<SdrfDriftFinding> Conflicts(
            Dictionary<string, Dictionary<string, Dictionary<string, List<string>>>> index,
            SdrfDriftKind kind)
        {
            foreach (var column in index)
            {
                foreach (var concept in column.Value.Where(c => c.Value.Count > 1))
                {
                    var variants = concept.Value
                        .Select(v => new SdrfDriftVariant(v.Key, v.Value.Distinct(StringComparer.Ordinal).Count(),
                            v.Value.Distinct(StringComparer.Ordinal).OrderBy(x => x, StringComparer.Ordinal).ToList()))
                        .OrderByDescending(v => v.Occurrences)
                        .ThenBy(v => v.Value, StringComparer.Ordinal)
                        .ToList();

                    yield return new SdrfDriftFinding(kind, concept.Key, column.Key, variants);
                }
            }
        }

        private static void Add(
            Dictionary<string, Dictionary<string, Dictionary<string, List<string>>>> index,
            string column, string concept, string variant, string label)
        {
            if (!index.TryGetValue(column, out var byConcept))
                index[column] = byConcept = new Dictionary<string, Dictionary<string, List<string>>>(StringComparer.Ordinal);
            if (!byConcept.TryGetValue(concept, out var byVariant))
                byConcept[concept] = byVariant = new Dictionary<string, List<string>>(StringComparer.Ordinal);
            if (!byVariant.TryGetValue(variant, out var documents))
                byVariant[variant] = documents = new List<string>();
            documents.Add(label);
        }

        private static void Track(Dictionary<string, HashSet<string>> index, string column, string label)
        {
            if (!index.TryGetValue(column, out var documents))
                index[column] = documents = new HashSet<string>(StringComparer.Ordinal);
            documents.Add(label);
        }

        private static bool IsReserved(string cell) =>
            SdrfValidator.ReservedWords.Any(w => string.Equals(cell, w, StringComparison.OrdinalIgnoreCase));

        /// <summary>
        /// Case- and whitespace-insensitive form, used to decide that two spellings are meant to be
        /// the same thing. Deliberately conservative: it does not touch punctuation or word order,
        /// so it will not claim that genuinely different values are variants of one another.
        /// </summary>
        private static string Normalize(string value) =>
            string.Join(' ', value.Split((char[])null, StringSplitOptions.RemoveEmptyEntries)).ToLowerInvariant();
    }
}
