using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using UglyToad.PdfPig.Content;
using UglyToad.PdfPig.Core;
using UglyToad.PdfPig.Fonts.Standard14Fonts;
using UglyToad.PdfPig.Writer;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Supplement tables as rows of strings (sdrf design SAMPLE-EVIDENCE.md, E2): .xlsx and .docx read as the zipped
    /// XML they are, delimited text read quote-aware, the real header row found under any title lines.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSupplementTableReader
    {
        private static string Dir => TestContext.CurrentContext.WorkDirectory;

        private static string Zip(string name, params (string Entry, string Xml)[] entries)
        {
            string path = Path.Combine(Dir, name);
            if (File.Exists(path)) File.Delete(path);
            using var zip = ZipFile.Open(path, ZipArchiveMode.Create);
            foreach (var (entry, xml) in entries)
            {
                using var w = new StreamWriter(zip.CreateEntry(entry).Open(), new UTF8Encoding(false));
                w.Write(xml);
            }
            return path;
        }

        private const string Ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
        private const string Rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";

        /// <summary>Two sheets: a title line over a sparse, mixed-type patient table; and a small second sheet.</summary>
        private static string Workbook() => Zip("patients.xlsx",
            ("xl/workbook.xml", $@"<workbook xmlns=""{Ns}"" xmlns:r=""{Rel}""><sheets>
                <sheet name=""Patients"" sheetId=""1"" r:id=""rId1""/><sheet name=""Batches"" sheetId=""2"" r:id=""rId2""/></sheets></workbook>"),
            ("xl/_rels/workbook.xml.rels", @"<Relationships xmlns=""http://schemas.openxmlformats.org/package/2006/relationships"">
                <Relationship Id=""rId1"" Target=""worksheets/sheet1.xml"" Type=""x""/><Relationship Id=""rId2"" Target=""/xl/worksheets/sheet2.xml"" Type=""x""/></Relationships>"),
            ("xl/sharedStrings.xml", $@"<sst xmlns=""{Ns}""><si><t>Table S1. Clinical data</t></si><si><t>ID</t></si><si><t>Gender (1F,0M)</t></si>
                <si><t>Age</t></si><si><r><t>Heart</t></r><r><t> failure</t></r></si><si><t>HumanHFpEF_1</t></si></sst>"),
            ("xl/worksheets/sheet1.xml", $@"<worksheet xmlns=""{Ns}""><sheetData>
                <row r=""1""><c r=""A1"" t=""s""><v>0</v></c></row>
                <row r=""2""><c r=""A2"" t=""s""><v>1</v></c><c r=""B2"" t=""s""><v>2</v></c><c r=""C2"" t=""s""><v>3</v></c><c r=""D2"" t=""s""><v>4</v></c></row>
                <row r=""3""><c r=""A3"" t=""s""><v>5</v></c><c r=""B3""><v>1</v></c><c r=""C3""><v>77</v></c><c r=""D3"" t=""b""><v>1</v></c></row>
                <row r=""5""><c r=""A5"" t=""inlineStr""><is><t>HumanControl_1</t></is></c><c r=""C5""><v>0.15</v></c></row>
                </sheetData></worksheet>"),
            ("xl/worksheets/sheet2.xml", $@"<worksheet xmlns=""{Ns}""><sheetData>
                <row r=""1""><c r=""A1"" t=""inlineStr""><is><t>batch</t></is></c><c r=""B1"" t=""inlineStr""><is><t>date</t></is></c></row>
                <row r=""2""><c r=""A2""><v>3.0</v></c><c r=""B2"" t=""str""><v>2024-01-02</v></c></row>
                </sheetData></worksheet>"));

        [Test]
        public void AWorkbookGivesEverySheetWithItsHeaderUnderTheTitle()
        {
            var tables = SupplementTableReader.Read(Workbook());

            Assert.That(tables.Select(t => t.Sheet), Is.EqualTo(new[] { "Patients", "Batches" }));
            var p = tables[0];
            Assert.That(p.Title, Is.EqualTo("Table S1. Clinical data"));
            Assert.That(p.Header, Is.EqualTo(new[] { "ID", "Gender (1F,0M)", "Age", "Heart failure" }), "rich text runs are joined");
            Assert.That(p.Rows, Has.Count.EqualTo(2), "the empty row 4 is dropped");
            Assert.That(p.Rows[0], Is.EqualTo(new[] { "HumanHFpEF_1", "1", "77", "TRUE" }));
            Assert.That(p.Rows[1], Is.EqualTo(new[] { "HumanControl_1", "", "0.15", "" }), "a sparse row is padded by column letter");
            Assert.That(tables[1].Rows[0], Is.EqualTo(new[] { "3", "2024-01-02" }), "3.0 is written as 3");
            Assert.That(p.Locator(0, 2), Is.EqualTo("patients.xlsx!Patients!R3C3"), "a locator names the sheet row and column");
        }

        [Test]
        public void ARowCapStopsALongSheet()
        {
            // The cap counts sheet rows, title and header included: rows 1-3 are read, row 5 is not.
            var tables = SupplementTableReader.Read(Workbook(), maxRows: 3);

            Assert.That(tables[0].Rows, Has.Count.EqualTo(1));
            Assert.That(tables[0].Truncated, Is.True);
        }

        [Test]
        public void AWordDocumentGivesItsTables()
        {
            const string w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main";
            string Cell(string t) => $"<w:tc><w:p><w:r><w:t>{t}</w:t></w:r></w:p></w:tc>";
            string path = Zip("s1.docx", ("word/document.xml", $@"<w:document xmlns:w=""{w}""><w:body>
                <w:p><w:r><w:t>Supplementary Table 2</w:t></w:r></w:p>
                <w:tbl><w:tr>{Cell("Sample")}{Cell("Sex")}</w:tr><w:tr>{Cell("C1")}{Cell("female")}</w:tr><w:tr>{Cell("C2")}{Cell("male")}</w:tr></w:tbl>
                <w:tbl><w:tr>{Cell("Antibody")}{Cell("Dilution")}</w:tr><w:tr>{Cell("anti-GFP")}{Cell("1:1000")}</w:tr></w:tbl>
                </w:body></w:document>"));

            var tables = SupplementTableReader.Read(path);

            Assert.That(tables.Select(t => t.Sheet), Is.EqualTo(new[] { "table 1", "table 2" }));
            Assert.That(tables[0].Header, Is.EqualTo(new[] { "Sample", "Sex" }));
            Assert.That(tables[0].Rows.Select(r => r[1]), Is.EqualTo(new[] { "female", "male" }));
        }

        [TestCase("meta.csv", "sample_id,group,note\nFL_c1,female,\"has, a comma\"\nFL_c2,male,\"says \"\"hi\"\"\"\n")]
        [TestCase("meta.tsv", "sample_id\tgroup\tnote\nFL_c1\tfemale\thas, a comma\nFL_c2\tmale\tsays \"hi\"\n")]
        [TestCase("meta.txt", "sample_id;group;note\nFL_c1;female;has, a comma\nFL_c2;male;says \"hi\"\n")]
        public void DelimitedTextIsReadQuoteAware(string name, string content)
        {
            string path = Path.Combine(Dir, name);
            File.WriteAllText(path, content);

            var t = SupplementTableReader.Read(path).Single();

            Assert.That(t.Header, Is.EqualTo(new[] { "sample_id", "group", "note" }));
            Assert.That(t.Rows[0][2], Is.EqualTo("has, a comma"));
            Assert.That(t.Rows[1][2], Is.EqualTo("says \"hi\""));
        }

        [Test]
        public void AnIsaTabFileReadsAsATable()
        {
            string path = Path.Combine(Dir, "s_study.txt");
            File.WriteAllText(path, "\"Source Name\"\t\"Characteristics[organism]\"\t\"Sample Name\"\n\"m1\"\t\"Mus musculus\"\t\"AJ_Heart_d00\"\n");

            var t = SupplementTableReader.Read(path).Single();

            Assert.That(t.Header, Is.EqualTo(new[] { "Source Name", "Characteristics[organism]", "Sample Name" }));
            Assert.That(t.Rows.Single(), Is.EqualTo(new[] { "m1", "Mus musculus", "AJ_Heart_d00" }));
        }

        [Test]
        public void AFormatWithNoReaderGivesNothingAndABrokenFileIsRefused()
        {
            string xls = Path.Combine(Dir, "legacy.xls");
            File.WriteAllText(xls, "not read");
            Assert.That(SupplementTableReader.Read(xls), Is.Empty);

            string pdf = Path.Combine(Dir, "broken.pdf");
            File.WriteAllText(pdf, "%PDF-1.4 and nothing else");
            Assert.Throws<MzLibException>(() => SupplementTableReader.Read(pdf));

            string broken = Path.Combine(Dir, "broken.xlsx");
            File.WriteAllText(broken, "not a zip");
            Assert.Throws<MzLibException>(() => SupplementTableReader.Read(broken));
        }

        /// <summary>A one-page PDF with text placed at (x, y) points, 10 pt Helvetica unless a size is given.</summary>
        private static string Pdf(string name, params (string Text, double X, double Y, double Size)[] words)
        {
            var builder = new PdfDocumentBuilder();
            var page = builder.AddPage(PageSize.A4);
            var font = builder.AddStandard14Font(Standard14Font.Helvetica);
            foreach (var (text, x, y, size) in words) page.AddText(text, size, new PdfPoint(x, y), font);
            string path = Path.Combine(Dir, name);
            File.WriteAllBytes(path, builder.Build());
            return path;
        }

        [Test]
        public void APdfTableIsReadFromWordPositionsWithItsCaptionAndAWrappedHeader()
        {
            // The shape of PXD021990's Table S7: a caption, a header cell wrapped over two lines ("Pediatric" above
            // "Age Group"), a superscript in a header (m2), and a cell wrapped onto the next line.
            string pdf = Pdf("s1.pdf",
                ("Supplementary Table S1. Patient characteristics", 50, 760, 10),
                ("Pediatric", 200, 730, 10),
                ("Patient", 50, 718, 10), ("Age", 130, 718, 10), ("Age Group", 200, 718, 10), ("Sex", 300, 718, 10), ("BSA (m", 360, 718, 10), ("2", 398, 722, 6), (")", 402, 718, 10),
                ("DB85", 50, 700, 10), ("9", 130, 700, 10), ("Child", 200, 700, 10), ("F", 300, 700, 10), ("0.95", 360, 700, 10),
                ("DB27", 50, 686, 10), ("17", 130, 686, 10), ("Adolescent", 200, 686, 10), ("M", 300, 686, 10), ("2.2", 360, 686, 10),
                ("(transplant)", 200, 674, 10),
                ("DB43", 50, 660, 10), ("7", 130, 660, 10), ("Infant", 200, 660, 10), ("M", 300, 660, 10), ("0.24", 360, 660, 10),
                ("This paragraph follows the table and is prose, set as one run of words across the page.", 50, 620, 10));

            var t = SupplementTableReader.Read(pdf).Single();

            Assert.That(t.Sheet, Is.EqualTo("page 1 table 1"));
            Assert.That(t.Title, Does.StartWith("Supplementary Table S1."));
            Assert.That(t.Header, Is.EqualTo(new[] { "Patient", "Age", "Pediatric Age Group", "Sex", "BSA (m 2 )" }));
            Assert.That(t.Rows.Select(r => r[0]), Is.EqualTo(new[] { "DB85", "DB27", "", "DB43" }), "the wrapped cell is its own row, as in the PDF");
            Assert.That(t.Rows[1], Is.EqualTo(new[] { "DB27", "17", "Adolescent", "M", "2.2" }));
            Assert.That(t.Rows[2][2], Is.EqualTo("(transplant)"));
            Assert.That(t.Locator(0, 2), Is.EqualTo("s1.pdf!page 1 table 1!R3C3"));
        }

        [Test]
        public void ProseInAPdfIsNotATable()
        {
            string pdf = Pdf("methods.pdf",
                ("Supplementary Methods", 50, 760, 12),
                ("Samples were collected from ten patients and ten donors,", 50, 740, 10),
                ("digested with trypsin and labelled with TMT 10-plex reagents.", 50, 726, 10),
                ("Each plex was fractionated into twelve fractions.", 50, 712, 10));

            Assert.That(SupplementTableReader.Read(pdf), Is.Empty);
        }
    }
}
