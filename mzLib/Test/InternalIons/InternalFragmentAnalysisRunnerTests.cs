using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Readers.InternalIons;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.InternalIons
{
    /// <summary>
    /// Integration harness for testing internal fragment ion extraction and analysis.
    /// Call RunAll() manually from any test or scratch runner.
    /// </summary>
    public static class InternalFragmentAnalysisRunnerTests
    {
        // Fill these in with your local test data paths
        private const string DefaultPsmTsvPath = @"FILL_IN";
        private const string DefaultRawFileFolder = @"FILL_IN";
        private const string DefaultOutputDirectory = @"FILL_IN";

        private static readonly char[] StandardAminoAcids =
            { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };

        private static readonly string[] SupportedExtensions = { ".raw", ".mzml", ".mgf" };

        public static void RunAll() => RunAll(DefaultPsmTsvPath, DefaultRawFileFolder, DefaultOutputDirectory);

        public static void RunAll(string psmTsvPath, string rawFileFolder, string outputDirectory)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("     INTERNAL FRAGMENT ION ANALYSIS HARNESS");
            Console.WriteLine("==================================================================");
            Console.WriteLine();
            Console.WriteLine($"PSM TSV:          {psmTsvPath}");
            Console.WriteLine($"Raw File Folder:  {rawFileFolder}");
            Console.WriteLine($"Output Directory: {outputDirectory}");
            Console.WriteLine();

            try
            {
                if (!File.Exists(psmTsvPath))
                    throw new FileNotFoundException($"PSM TSV file not found: {psmTsvPath}");
                if (!Directory.Exists(rawFileFolder))
                    throw new DirectoryNotFoundException($"Raw file folder not found: {rawFileFolder}");
                if (!Directory.Exists(outputDirectory))
                    Directory.CreateDirectory(outputDirectory);

                var psms = Step1_LoadPsms(psmTsvPath);
                var msDataFiles = Step2_LoadRawFiles(rawFileFolder, psms);
                var internalIons = Step3_ExtractInternalFragments(psms, msDataFiles);
                string outputPath = Step4_WriteOutputTsv(internalIons, outputDirectory);
                Step5_RoundTripValidation(internalIons, outputPath);
                Step6_AminoAcidTerminusEnrichment(internalIons);

                Console.WriteLine();
                Console.WriteLine("==================================================================");
                Console.WriteLine("ALL STEPS COMPLETED SUCCESSFULLY");
                Console.WriteLine("==================================================================");
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                Console.WriteLine("==================================================================");
                Console.WriteLine("FAILURE");
                Console.WriteLine("==================================================================");
                Console.WriteLine($"Exception: {ex.GetType().FullName}");
                Console.WriteLine($"Message: {ex.Message}");
                Console.WriteLine($"Stack Trace:\n{ex.StackTrace}");
                throw;
            }
        }

        private static List<PsmFromTsv> Step1_LoadPsms(string psmTsvPath)
        {
            Console.WriteLine("=== STEP 1: Load PSMs ===");
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var warnings);
            Console.WriteLine($"Total PSMs loaded: {psms.Count:N0}");

            var sourceFiles = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToList();
            Console.WriteLine($"Unique source files: {sourceFiles.Count}");

            var psmsWithInternal = psms
                .Where(p => p.MatchedIons != null && p.MatchedIons.Any(ion => IsInternalAnnotation(ion.Annotation)))
                .ToList();
            Console.WriteLine($"PSMs with internal fragments: {psmsWithInternal.Count:N0}");

            foreach (var psm in psmsWithInternal.Take(3))
            {
                var annotations = psm.MatchedIons
                    .Where(ion => IsInternalAnnotation(ion.Annotation))
                    .Select(ion => ion.Annotation)
                    .Take(3);
                Console.WriteLine($"  Scan {psm.Ms2ScanNumber}: {psm.BaseSeq} -> {string.Join(", ", annotations)}");
            }
            Console.WriteLine();
            return psms;
        }

        private static bool IsInternalAnnotation(string annotation)
        {
            if (string.IsNullOrEmpty(annotation)) return false;
            if (annotation.Contains("int", StringComparison.OrdinalIgnoreCase)) return true;
            int start = annotation.IndexOf('[');
            int end = annotation.IndexOf(']');
            return start >= 0 && end > start && annotation.Substring(start, end - start).Contains('-');
        }

        private static Dictionary<string, MsDataFile> Step2_LoadRawFiles(string rawFileFolder, List<PsmFromTsv> psms)
        {
            Console.WriteLine("=== STEP 2: Load Raw Files ===");
            var rawFiles = SupportedExtensions
                .SelectMany(ext => Directory.GetFiles(rawFileFolder, $"*{ext}", SearchOption.TopDirectoryOnly))
                .ToList();
            Console.WriteLine($"Found {rawFiles.Count} raw files");

            var neededFiles = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToHashSet(StringComparer.OrdinalIgnoreCase);
            var msDataFiles = new Dictionary<string, MsDataFile>(StringComparer.OrdinalIgnoreCase);
            int totalMs2 = 0;

            foreach (var rawFile in rawFiles)
            {
                string name = Path.GetFileNameWithoutExtension(rawFile);
                if (!neededFiles.Contains(name)) continue;

                Console.WriteLine($"  Loading {Path.GetFileName(rawFile)}...");
                try
                {
                    string ext = Path.GetExtension(rawFile).ToLowerInvariant();
                    MsDataFile file = ext switch
                    {
                        ".mzml" => Mzml.LoadAllStaticData(rawFile),
                        ".raw" => ThermoRawFileReader.LoadAllStaticData(rawFile),
                        ".mgf" => Mgf.LoadAllStaticData(rawFile),
                        _ => throw new NotSupportedException(ext)
                    };
                    msDataFiles[name] = file;
                    int ms2 = file.GetAllScansList().Count(s => s.MsnOrder == 2);
                    totalMs2 += ms2;
                    Console.WriteLine($"    {file.GetAllScansList().Count} scans ({ms2} MS2)");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"    ERROR: {ex.Message}");
                }
            }
            Console.WriteLine($"Total MS2 scans: {totalMs2:N0}");
            Console.WriteLine();
            return msDataFiles;
        }

        private static List<InternalFragmentIon> Step3_ExtractInternalFragments(
            List<PsmFromTsv> psms, Dictionary<string, MsDataFile> msDataFiles)
        {
            Console.WriteLine("=== STEP 3: Extract Internal Fragment Features ===");
            var allIons = new List<InternalFragmentIon>();

            foreach (var group in psms.GroupBy(p => p.FileNameWithoutExtension))
            {
                if (!msDataFiles.TryGetValue(group.Key, out var file)) continue;
                var ions = InternalFragmentFeatureExtractor.ExtractFromPsms(group.ToList(), file);
                allIons.AddRange(ions);
                if (ions.Count > 0)
                    Console.WriteLine($"  {group.Key}: {ions.Count:N0} internal ions");
            }

            Console.WriteLine($"Total internal ions: {allIons.Count:N0}");

            if (allIons.Count > 0)
            {
                int proline = allIons.Count(i => i.HasProlineAtEitherTerminus);
                int aspartate = allIons.Count(i => i.HasAspartateAtEitherTerminus);
                Console.WriteLine($"With Proline terminus: {proline:N0} ({100.0 * proline / allIons.Count:F1}%)");
                Console.WriteLine($"With Aspartate terminus: {aspartate:N0} ({100.0 * aspartate / allIons.Count:F1}%)");

                var errors = allIons.Select(i => i.MassErrorPpm).Where(e => !double.IsNaN(e)).ToList();
                if (errors.Count > 0)
                {
                    double mean = errors.Average();
                    double std = Math.Sqrt(errors.Select(e => Math.Pow(e - mean, 2)).Average());
                    Console.WriteLine($"Mass Error (ppm): Mean={mean:F4}, StdDev={std:F4}");
                }
            }
            Console.WriteLine();
            return allIons;
        }

        private static string Step4_WriteOutputTsv(List<InternalFragmentIon> ions, string outputDirectory)
        {
            Console.WriteLine("=== STEP 4: Write Output TSV ===");
            string outputPath = Path.Combine(outputDirectory, "InternalFragmentIons.tsv");
            InternalFragmentTsvWriter.WriteToTsv(ions, outputPath);
            Console.WriteLine($"Output: {outputPath}");
            Console.WriteLine($"Rows: {ions.Count:N0}");
            Console.WriteLine();
            return outputPath;
        }

        private static void Step5_RoundTripValidation(List<InternalFragmentIon> original, string outputPath)
        {
            Console.WriteLine("=== STEP 5: Round-trip Validation ===");
            if (original.Count == 0)
            {
                Console.WriteLine("Skipping (no ions)");
                Console.WriteLine();
                return;
            }

            var readBack = InternalFragmentTsvWriter.ReadFromTsv(outputPath);
            Console.WriteLine($"Original: {original.Count}, Read back: {readBack.Count}");

            if (original.Count != readBack.Count)
                throw new InvalidOperationException($"Count mismatch: {original.Count} vs {readBack.Count}");

            var o = original[0];
            var r = readBack[0];
            if (o.PeptideSequence != r.PeptideSequence || o.InternalSequence != r.InternalSequence)
                throw new InvalidOperationException("Sequence mismatch on first ion");
            if (Math.Abs(o.NormalizedIntensity - r.NormalizedIntensity) > 1e-6)
                throw new InvalidOperationException("NormalizedIntensity mismatch");

            Console.WriteLine("Round-trip validation PASSED");
            Console.WriteLine();
        }

        private static void Step6_AminoAcidTerminusEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 6: Amino Acid Terminus Enrichment ===");
            if (ions.Count == 0)
            {
                Console.WriteLine("Skipping (no ions)");
                return;
            }

            string allSeq = string.Concat(ions.Select(i => i.InternalSequence));
            var bgFreq = StandardAminoAcids.ToDictionary(aa => aa, aa => (double)allSeq.Count(c => c == aa) / allSeq.Length);
            var termCount = StandardAminoAcids.ToDictionary(aa => aa, _ => 0);

            foreach (var ion in ions)
            {
                if (string.IsNullOrEmpty(ion.InternalSequence)) continue;
                char n = ion.InternalSequence[0];
                char c = ion.InternalSequence[^1];
                if (termCount.ContainsKey(n)) termCount[n]++;
                if (termCount.ContainsKey(c)) termCount[c]++;
            }

            int totalPos = ions.Count * 2;
            Console.WriteLine("+----+----------+----------+------------+-------+");
            Console.WriteLine("| AA | Observed | Expected | Enrichment | Flag  |");
            Console.WriteLine("+----+----------+----------+------------+-------+");

            var data = StandardAminoAcids
                .Select(aa => (aa, obs: termCount[aa], exp: bgFreq[aa] * totalPos, enrich: bgFreq[aa] > 0 ? termCount[aa] / (bgFreq[aa] * totalPos) : 0))
                .OrderByDescending(x => x.enrich);

            foreach (var (aa, obs, exp, enrich) in data)
            {
                string flag = enrich > 1.5 ? " ***" : "";
                Console.WriteLine($"|  {aa} | {obs,8} | {exp,8:F1} | {enrich,10:F3} |{flag,-6}|");
            }
            Console.WriteLine("+----+----------+----------+------------+-------+");
            Console.WriteLine("*** = Enrichment > 1.5");
            Console.WriteLine();
        }
    }

    [TestFixture]
    public class InternalFragmentTests
    {
        [Test]
        [Explicit("Manual integration harness - requires local data files")]
        public void RunInternalFragmentAnalysis()
        {
            InternalFragmentAnalysisRunnerTests.RunAll(
                psmTsvPath: @"F:\JurakatMultiprotease\2026-02-17-16-40-42\Task1-SearchTask\AllPeptides.psmtsv",
                rawFileFolder: @"F:\JurakatMultiprotease",
                outputDirectory: @"E:\Projects\internalIons\mzLibReportsForClaude"
            );
        }
    }
}