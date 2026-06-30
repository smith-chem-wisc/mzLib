using System.Collections.Generic;
using System.IO;
using NUnit.Framework;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>One annotated manuscript example. Columns mirror manuscript-examples.tsv.</summary>
    internal record CorpusRow(string Id, string ProformaString, string Source, string Section,
        string ComplianceLevel, string Feature, string CanonicalForm, bool Valid, string Notes);

    /// <summary>
    /// Loads the manuscript-example corpus (the authoritative ProForma test set) once for all
    /// ProForma test fixtures. Source enumeration: corpus/EXAMPLE-INDEX.md.
    /// </summary>
    internal static class ProFormaTestCorpus
    {
        public static string TsvPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ProForma\manuscript-examples.tsv");

        public static List<CorpusRow> Load()
        {
            var rows = new List<CorpusRow>();
            var lines = File.ReadAllLines(TsvPath);
            for (int i = 1; i < lines.Length; i++) // skip header
            {
                if (lines[i].Length == 0) continue;
                var c = lines[i].Split('\t');
                // Throw a plain data error (not an NUnit Assert): Load() runs from [TestCaseSource] at
                // discovery time, where an Assert would surface as an unreadable fixture-load failure.
                if (c.Length != 9)
                    throw new InvalidDataException(
                        $"Row {i + 1} of {TsvPath} does not have 9 tab-separated columns (has {c.Length}).");
                rows.Add(new CorpusRow(c[0], c[1], c[2], c[3], c[4], c[5], c[6], bool.Parse(c[7]), c[8]));
            }
            return rows;
        }
    }
}
