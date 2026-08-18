using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for pooling SDRF documents and for detecting inconsistency between them.
    ///
    /// The distinction being tested is the point of the whole component: every document here is
    /// individually VALID -- SdrfValidator passes all of them -- and the pooled result is still
    /// unusable, because validity is a property of one file and comparability is a property of the
    /// relationship between files.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfCollection
    {
        private static string FixtureDir => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes");

        private static SdrfDocument Doc(string[] columns, params string[][] rows)
        {
            var header = new SdrfHeader(columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        // ---------------- merge ----------------

        [Test]
        public void Merge_UnionsColumnsAndRecordsProvenance()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "characteristics[organism]" },
                new[] { "S1", "run 1", "t", "homo sapiens" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]" },
                new[] { "S1", "run 1", "t", "NT=Q Exactive;AC=MS:1001911" });

            var merged = new SdrfCollection(new[] { a, b }, new[] { "PXD000001", "PXD000002" }).Merge();

            Assert.That(merged.Header.Contains("characteristics[organism]"), Is.True);
            Assert.That(merged.Header.Contains("comment[instrument]"), Is.True);
            Assert.That(merged.Header.Contains(SdrfCollection.SourceDocumentColumn), Is.True);
            Assert.That(merged.Results.Count, Is.EqualTo(2));

            Assert.That(merged.Results[0][SdrfCollection.SourceDocumentColumn], Is.EqualTo("PXD000001"));
            Assert.That(merged.Results[1][SdrfCollection.SourceDocumentColumn], Is.EqualTo("PXD000002"));

            // A column the other document lacks says so explicitly rather than being blank.
            Assert.That(merged.Results[1]["characteristics[organism]"], Is.EqualTo("not available"));
            Assert.That(merged.Results[0]["comment[instrument]"], Is.EqualTo("not available"));
        }

        [Test]
        public void Merge_KeepsSdrfBlockOrder()
        {
            var a = Doc(new[] { "comment[instrument]", "factor value[disease]", "source name", "assay name", "technology type", "characteristics[organism]" },
                new[] { "i", "d", "S1", "run 1", "t", "o" });

            var merged = new SdrfCollection(new[] { a }, new[] { "X" }).Merge();
            var names = merged.Header.ToList();

            Assert.That(names.IndexOf("source name"), Is.LessThan(names.IndexOf("characteristics[organism]")));
            Assert.That(names.IndexOf("characteristics[organism]"), Is.LessThan(names.IndexOf("assay name")));
            Assert.That(names.IndexOf("assay name"), Is.LessThan(names.IndexOf("comment[instrument]")));
            Assert.That(names.IndexOf("comment[instrument]"), Is.LessThan(names.IndexOf("factor value[disease]")));
        }

        [Test]
        public void Merge_PreservesRepeatedColumnMultiplicity()
        {
            // comment[modification parameters] repeats up to 8 times in the corpus. A merge that
            // deduplicated column names would silently drop all but the first modification.
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[modification parameters]", "comment[modification parameters]" },
                new[] { "S1", "run 1", "t", "NT=Oxidation;AC=UNIMOD:35", "NT=Acetyl;AC=UNIMOD:1" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[modification parameters]" },
                new[] { "S2", "run 2", "t", "NT=Phospho;AC=UNIMOD:21" });

            var merged = new SdrfCollection(new[] { a, b }, new[] { "A", "B" }).Merge();

            Assert.That(merged.Header.IndexesOf("comment[modification parameters]").Count, Is.EqualTo(2));
            Assert.That(merged.Results[0].All("comment[modification parameters]"),
                Is.EqualTo(new[] { "NT=Oxidation;AC=UNIMOD:35", "NT=Acetyl;AC=UNIMOD:1" }));
            Assert.That(merged.Results[1].All("comment[modification parameters]"),
                Is.EqualTo(new[] { "NT=Phospho;AC=UNIMOD:21", "not available" }));
        }

        [Test]
        public void Merge_IsWritableAsSdrf()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type" }, new[] { "S1", "run 1", "t" });
            var merged = new SdrfCollection(new[] { a }, new[] { "A" }).Merge();

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"merged_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                merged.WriteResults(outputPath);
                var reloaded = new SdrfDocument(outputPath);
                Assert.That(reloaded.Header.Contains(SdrfCollection.SourceDocumentColumn), Is.True);
                Assert.That(reloaded.Results.Count, Is.EqualTo(1));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        // ---------------- drift ----------------

        [Test]
        public void Drift_SameNameDifferentAccession_IsDetected()
        {
            // The serious one, and a live problem: mzLib's proteases.tsv gives plain trypsin
            // MS:1001313 (Trypsin/P) while PSI-MS plain Trypsin is MS:1001251. Two documents
            // annotating "Trypsin" from different sources split silently.
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[cleavage agent details]" },
                new[] { "S1", "run 1", "t", "NT=Trypsin;AC=MS:1001251" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[cleavage agent details]" },
                new[] { "S2", "run 2", "t", "NT=Trypsin;AC=MS:1001313" });

            var findings = SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" }));

            var finding = findings.FirstOrDefault(f => f.Kind == SdrfDriftKind.NameAccessionConflict);
            Assert.That(finding, Is.Not.Null);
            Assert.That(finding!.Concept, Is.EqualTo("Trypsin"));
            Assert.That(finding.Variants.Select(v => v.Value),
                Is.EquivalentTo(new[] { "MS:1001251", "MS:1001313" }));
        }

        [Test]
        public void Drift_SameAccessionDifferentName_IsDetected()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]" },
                new[] { "S1", "run 1", "t", "NT=Q Exactive;AC=MS:1001911" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]" },
                new[] { "S2", "run 2", "t", "NT=q exactive;AC=MS:1001911" });

            var findings = SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" }));

            var finding = findings.FirstOrDefault(f => f.Kind == SdrfDriftKind.AccessionNameConflict);
            Assert.That(finding, Is.Not.Null);
            Assert.That(finding!.Concept, Is.EqualTo("MS:1001911"));
        }

        [Test]
        public void Drift_TermVersusFreeText_IsDetected()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "characteristics[organism]" },
                new[] { "S1", "run 1", "t", "NT=Homo sapiens;AC=NCBITaxon:9606" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "characteristics[organism]" },
                new[] { "S2", "run 2", "t", "Homo sapiens" });

            var findings = SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" }));

            Assert.That(findings.Any(f => f.Kind == SdrfDriftKind.MixedTermAndFreeText
                                          && f.Column == "characteristics[organism]"), Is.True);
        }

        [Test]
        public void Drift_FreeTextCaseVariant_IsDetected()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[label]" },
                new[] { "S1", "run 1", "t", "label free sample" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[label]" },
                new[] { "S2", "run 2", "t", "Label Free Sample" });

            var findings = SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" }));

            var finding = findings.FirstOrDefault(f => f.Kind == SdrfDriftKind.ValueCaseVariant);
            Assert.That(finding, Is.Not.Null);
            Assert.That(finding!.Majority.Occurrences, Is.EqualTo(1));
        }

        [Test]
        public void Drift_ColumnNameVariant_IsDetected()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "characteristics[organism part]" },
                new[] { "S1", "run 1", "t", "liver" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "Characteristics[Organism Part]" },
                new[] { "S2", "run 2", "t", "liver" });

            var findings = SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" }));

            Assert.That(findings.Any(f => f.Kind == SdrfDriftKind.ColumnNameVariant), Is.True);
        }

        [Test]
        public void Drift_ConsistentDocuments_ProduceNoFindings()
        {
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]" },
                new[] { "S1", "run 1", "t", "NT=Q Exactive;AC=MS:1001911" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]" },
                new[] { "S2", "run 2", "t", "NT=Q Exactive;AC=MS:1001911" });

            Assert.That(SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" })), Is.Empty);
        }

        [Test]
        public void Drift_PerRowColumnsAreIgnored()
        {
            // source name and comment[data file] are unique per row by design. Reporting them would
            // bury every real finding.
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[data file]" },
                new[] { "Sample 1", "run 1", "t", "a.raw" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[data file]" },
                new[] { "sample 1", "run 1", "t", "A.RAW" });

            Assert.That(SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" })), Is.Empty);
        }

        [Test]
        public void Drift_ValidDocumentsCanStillBeIncomparable()
        {
            // The whole thesis, stated as a test: per-file validity says nothing about whether two
            // files can be pooled. Both documents therefore carry every REQUIRED column, so each is
            // genuinely valid on its own -- the assertion below would be vacuous otherwise.
            string[] names =
            {
                "source name", "characteristics[organism]", "characteristics[organism part]",
                "characteristics[biological replicate]", "assay name", "technology type",
                "comment[proteomics data acquisition method]", "comment[label]", "comment[instrument]",
                "comment[cleavage agent details]", "comment[fraction identifier]",
                "comment[technical replicate]", "comment[data file]"
            };
            var a = Doc(names, new[]
            {
                "S1", "homo sapiens", "liver", "1", "run 1", "t",
                "NT=Data-dependent acquisition;AC=PRIDE:0000627", "label free sample",
                "NT=Q Exactive;AC=MS:1001911", "NT=Trypsin;AC=MS:1001251", "1", "1", "a.raw"
            });
            var b = Doc(names, new[]
            {
                "S2", "homo sapiens", "liver", "1", "run 2", "t",
                "NT=Data-dependent acquisition;AC=PRIDE:0000627", "label free sample",
                "NT=Q Exactive;AC=MS:1001911", "NT=Trypsin;AC=MS:1001313", "1", "1", "b.raw"
            });

            Assert.That(SdrfValidator.Validate(a).IsValid, Is.True);
            Assert.That(SdrfValidator.Validate(b).IsValid, Is.True);
            Assert.That(SdrfDriftLint.Analyze(new SdrfCollection(new[] { a, b }, new[] { "A", "B" })), Is.Not.Empty);
        }

        // ---------------- cell formatting ----------------

        [Test]
        public void ToCell_RoundTripsThroughTryParseTerm()
        {
            var term = new MzLibUtil.CvParam("UNIMOD", "UNIMOD:35", "Oxidation", "");
            string cell = SdrfCell.ToCell(term, ("TA", "M"), ("MT", "Variable"));

            Assert.That(cell, Is.EqualTo("NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable"));
            Assert.That(SdrfCell.IsTerm(cell), Is.True);
            Assert.That(SdrfCell.TryParseTerm(cell, out var parsed), Is.True);
            Assert.That(parsed!.Accession, Is.EqualTo("UNIMOD:35"));
            Assert.That(parsed.Name, Is.EqualTo("Oxidation"));

            var pairs = SdrfCell.ParseKeyValues(cell);
            Assert.That(pairs["TA"], Is.EqualTo("M"));
            Assert.That(pairs["MT"], Is.EqualTo("Variable"));
        }

        [Test]
        public void TryParseTerm_HandlesTheHalfTermsTheCorpusActuallyWrites()
        {
            // A name with no accession is still a term, and the CV label must stay empty rather
            // than being guessed, because there is no prefix to derive one from.
            Assert.That(SdrfCell.TryParseTerm("NT=Orbitrap", out var named), Is.True);
            Assert.That(named!.Name, Is.EqualTo("Orbitrap"));
            Assert.That(named.Accession, Is.Empty);
            Assert.That(named.CvLabel, Is.Empty);

            // The mirror case: an accession with no name is still a term, and the label comes off
            // the prefix.
            Assert.That(SdrfCell.TryParseTerm("AC=MS:1001911", out var accessionOnly), Is.True);
            Assert.That(accessionOnly!.Name, Is.Empty);
            Assert.That(accessionOnly.CvLabel, Is.EqualTo("MS"));
        }

        [Test]
        public void BareN_IsFreeTextNotAName()
        {
            // "N" is not a key in the SDRF grammar -- the specification defines NT as the only name
            // key. 45 corpus cells write "N=Orbitrap" in comment[ms2 analyzer type] (27 in
            // PXD039582, 18 in PXD039585), and an earlier revision accepted them as terms. That
            // decoded one submitter's typo AS the grammar, which also hid it from the drift lint,
            // whose whole job is to report exactly this.
            Assert.That(SdrfCell.IsTerm("N=Orbitrap"), Is.False);
            Assert.That(SdrfCell.TryParseTerm("N=Orbitrap", out var term), Is.False);
            Assert.That(term, Is.Null);
            Assert.That(SdrfCell.ParseKeyValues("N=Orbitrap"), Is.Empty);

            // And it is not silently readmitted after a legitimate leading key either.
            Assert.That(SdrfCell.TryParseTerm("AC=MS:1000484;N=Orbitrap", out var partial), Is.True);
            Assert.That(partial!.Name, Is.Empty, "N must not be promoted to the name");
        }

        [Test]
        public void SnLeadsPooledSampleCellsAndStaysInTheGrammar()
        {
            // SN is easy to mistake for drift because it is absent from the modification-style key
            // list, but the specification defines it for characteristics[pooled sample]
            // ("SN=sample1;SN=sample2", README.adoc:352 and :362). All 383 leading SN cells in the
            // corpus are of exactly that shape, so it stays -- unlike N, which the spec never
            // defines.
            Assert.That(SdrfCell.IsTerm("SN=sample_1;SN=sample_2"), Is.True);
        }

        [Test]
        public void LeadingKeysAreASubsetOfKnownKeys()
        {
            // KnownKeys is the source of truth; a key cannot enter the grammar through the
            // leading-key set without being justified there first. Had this existed, it would have
            // caught "N" on its own -- it was in both sets and in neither the specification nor any
            // deliberate decision.
            Assert.That(SdrfCell.KnownLeadingKeys, Is.SubsetOf(SdrfCell.KnownKeys));
        }

        [Test]
        public void TryParseTerm_RefusesToPromoteAKeyValueCellThatNamesNoTerm()
        {
            // characteristics[pooled sample] carries up to 45 repeated SN= keys. ParseKeyValues is
            // last-wins, so promoting these produced a Name that was one arbitrary member of the
            // list with an empty accession -- and because callers branch on this method, 538 cells
            // left the free-text index without entering the accession index, going invisible to
            // every kind of drift analysis at once. Neither NT nor AC means not a term.
            const string pooled = "SN=OSL.53E;SN=OSL.567";

            Assert.That(SdrfCell.TryParseTerm(pooled, out var term), Is.False);
            Assert.That(term, Is.Null);
            Assert.That(SdrfCell.IsTerm(pooled), Is.True,
                "IsTerm still reports the key=value grammar; only TryParseTerm declines to promote it");
        }

        [Test]
        public void ToCell_WritesTheAccessionVerbatim()
        {
            // Deliberately does NOT upper-case or add a missing prefix. The corpus is full of that
            // drift (bare "4", "Unimod:35", bare "1001251") and laundering it here would hide the
            // bug and make the drift lint's job impossible.
            var sloppy = new MzLibUtil.CvParam("", "4", "Carbamidomethyl", "");
            Assert.That(SdrfCell.ToCell(sloppy), Is.EqualTo("NT=Carbamidomethyl;AC=4"));
        }

        [Test]
        public void ToCell_RejectsASeparatorInsideAValue()
        {
            var term = new MzLibUtil.CvParam("MS", "MS:1001911", "Q\tExactive", "");
            Assert.That(() => SdrfCell.ToCell(term), Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void ToCell_RequiresANameOrAnAccession()
        {
            Assert.That(() => SdrfCell.ToCell(new MzLibUtil.CvParam("", "", "", "")),
                Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void Drift_FindingOrderIsDeterministic()
        {
            // Findings are built by walking Dictionary instances, which guarantee no enumeration
            // order. Unordered output is harmless in a console report and fatal the moment the
            // report is diffed between runs to see whether the corpus is drifting.
            var a = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]", "comment[label]" },
                new[] { "S1", "run 1", "t", "NT=Q Exactive;AC=MS:1001911", "label free sample" });
            var b = Doc(new[] { "source name", "assay name", "technology type", "comment[instrument]", "comment[label]" },
                new[] { "S2", "run 2", "t", "NT=q exactive;AC=MS:1001911", "Label Free Sample" });

            var collection = new SdrfCollection(new[] { a, b }, new[] { "A", "B" });

            var first = SdrfDriftLint.Analyze(collection).Select(f => f.ToString()).ToList();
            for (int i = 0; i < 5; i++)
            {
                var again = SdrfDriftLint.Analyze(collection).Select(f => f.ToString()).ToList();
                Assert.That(again, Is.EqualTo(first));
            }
            Assert.That(first, Is.Not.Empty);
        }

        [Test]
        public void Merge_DoesNotDuplicateAnExistingProvenanceColumn()
        {
            // Merging a merged document is the obvious incremental workflow. Appending the
            // provenance column unconditionally duplicated it, overwrote the original document's
            // provenance, and left the new one "not available" on every row.
            var a = Doc(new[] { "source name", "assay name", "technology type" }, new[] { "S1", "run 1", "t" });
            var b = Doc(new[] { "source name", "assay name", "technology type" }, new[] { "S2", "run 2", "t" });

            var once = new SdrfCollection(new[] { a, b }, new[] { "A", "B" }).Merge();
            var twice = new SdrfCollection(new[] { once }, new[] { "MERGED" }).Merge();

            Assert.That(twice.Header.IndexesOf(SdrfCollection.SourceDocumentColumn).Count, Is.EqualTo(1));
            Assert.That(twice.Results.Select(r => r[SdrfCollection.SourceDocumentColumn]),
                Has.None.EqualTo("not available"));
        }

        [Test]
        public void Collection_FromFiles_LabelsByAccession()
        {
            var collection = SdrfCollection.FromFiles(new[]
            {
                Path.Combine(FixtureDir, "PXD000070.sdrf.tsv"),
                Path.Combine(FixtureDir, "PXD026824.sdrf.tsv")
            });

            Assert.That(collection.Count, Is.EqualTo(2));

            // Labels are folder-qualified, not bare file names: every search writes its SDRF into
            // its own output folder, so N re-searches of one experiment share a base name and bare
            // names would collide.
            Assert.That(collection.Labels.Count, Is.EqualTo(2));
            Assert.That(collection.Labels[0], Does.EndWith("PXD000070"));
            Assert.That(collection.Labels[1], Does.EndWith("PXD026824"));
            Assert.That(collection.Labels[0], Is.Not.EqualTo(collection.Labels[1]));
            Assert.That(collection.Merge().Results.Count,
                Is.EqualTo(collection[0].Results.Count + collection[1].Results.Count));
        }

        /// <summary>
        /// Reports the drift the community's own curated corpus contains. [Explicit]; needs
        /// MZLIB_SDRF_CORPUS.
        ///
        /// Not a pass/fail check on anyone's annotations -- it is the vocabulary survey that tells
        /// us which spelling to adopt so our own documents stay joinable to public data. The
        /// majority variant of each finding is the answer.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void CorpusDriftReport()
        {
            string corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories);
            var collection = SdrfCollection.FromFiles(files);
            var findings = SdrfDriftLint.Analyze(collection);

            TestContext.Progress.WriteLine($"documents : {collection.Count}");
            TestContext.Progress.WriteLine($"findings  : {findings.Count}");
            foreach (var group in findings.GroupBy(f => f.Kind).OrderByDescending(g => g.Count()))
                TestContext.Progress.WriteLine($"  {group.Key,-24} {group.Count()}");

            // Ranking by variant count is dominated by a handful of pathological files -- e.g.
            // PXD029009-DIA gives the same instrument a different accession on every row. Ranking by
            // how many documents the MINORITY spellings touch surfaces the drift that actually
            // costs us joins, which is the guidance we want.
            TestContext.Progress.WriteLine("");
            TestContext.Progress.WriteLine("most impactful (documents using a non-majority spelling):");
            foreach (var finding in findings
                         .Select(f => (f, Minority: f.Variants.Skip(1).Sum(v => v.Occurrences)))
                         .Where(x => x.Minority > 0)
                         .OrderByDescending(x => x.Minority)
                         .Take(20))
            {
                var top = finding.f.Variants.Take(4).Select(v => $"{v.Value} x{v.Occurrences}");
                TestContext.Progress.WriteLine(
                    $"  [{finding.Minority,4}] {finding.f.Kind} '{finding.f.Concept}'" +
                    (finding.f.Column is null ? "" : $" in [{finding.f.Column}]") +
                    " :: " + string.Join(" | ", top));
            }

            TestContext.Progress.WriteLine("");
            TestContext.Progress.WriteLine("worst offenders (most competing spellings):");
            foreach (var finding in findings.OrderByDescending(f => f.Variants.Count).Take(8))
                TestContext.Progress.WriteLine("  " + Truncate(finding.ToString(), 200));

            Assert.That(findings, Is.Not.Null);
        }

        private static string Truncate(string value, int max) =>
            value.Length <= max ? value : value.Substring(0, max) + " ...";
    }
}
