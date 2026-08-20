using System.Collections;

namespace Readers
{
    /// <summary>
    /// The ordered column names of an SDRF-Proteomics document.
    ///
    /// This is a list, not a dictionary, and that is deliberate. SDRF column names are DATA, not a
    /// fixed schema, and two properties of real files rule a dictionary out:
    ///
    ///   1. Names repeat. In the 1,236-file curated corpus at bigbio/sdrf-annotated-datasets,
    ///      649 files repeat "comment[modification parameters]" (up to 8 times in one file),
    ///      258 repeat "comment[sdrf template]", and 102 repeat "comment[cleavage agent details]".
    ///      One file repeats an EMPTY name 23 times (trailing tabs).
    ///   2. Order is semantically meaningful. The specification fixes the block order
    ///      (sample metadata -> data file metadata -> factor values), so reordering columns
    ///      changes the document.
    ///
    /// Names are also stored verbatim, never case-normalized. The specification says column names
    /// are lowercase and case-sensitive, but the corpus contains "comment[MS min charge]",
    /// "Material Type", and "Factor Value[organism part]". Normalizing here would silently rewrite
    /// a user's file; reporting the deviation is the validator's job, not the reader's.
    /// </summary>
    public sealed class SdrfHeader : IReadOnlyList<string>
    {
        private readonly List<string> _names;

        public SdrfHeader(IEnumerable<string> names)
        {
            _names = names?.ToList() ?? throw new ArgumentNullException(nameof(names));
        }

        public string this[int index] => _names[index];

        public int Count => _names.Count;

        /// <summary>
        /// Every column index carrying this name, in document order. Empty when the column is absent.
        /// Callers that want a single value should use <see cref="IndexOf"/>; callers dealing with a
        /// multi-cardinality column (modification parameters, cleavage agent details) want all of them.
        /// Comparison is ordinal and case-SENSITIVE, matching the specification.
        /// </summary>
        public IReadOnlyList<int> IndexesOf(string columnName)
        {
            var found = new List<int>();
            for (int i = 0; i < _names.Count; i++)
            {
                if (string.Equals(_names[i], columnName, StringComparison.Ordinal))
                    found.Add(i);
            }
            return found;
        }

        /// <summary>
        /// The index of the first column with this name, or -1 if absent.
        /// </summary>
        public int IndexOf(string columnName) => _names.FindIndex(n => string.Equals(n, columnName, StringComparison.Ordinal));

        public bool Contains(string columnName) => IndexOf(columnName) >= 0;

        public IEnumerator<string> GetEnumerator() => _names.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        public override string ToString() => string.Join('\t', _names);
    }
}
