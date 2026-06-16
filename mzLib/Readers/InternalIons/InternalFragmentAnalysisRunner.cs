using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Readers.InternalIons
{
    /// <summary>
    /// Orchestrates the full internal fragment ion analysis workflow.
    /// </summary>
    public static class InternalFragmentAnalysisRunner
    {
        /// <summary>
        /// Runs the complete internal fragment extraction pipeline.
        /// </summary>
        /// <param name="psmTsvPath">Path to MetaMorpheus PSM TSV file.</param>
        /// <param name="rawFileFolder">Folder containing raw spectra files (.raw, .mzML).</param>
        /// <param name="outputDirectory">Directory for output files.</param>
        /// <param name="defaultCollisionEnergy">Default collision energy to use if scan metadata is unavailable.</param>
        /// <returns>List of extracted internal fragment ions.</returns>
        public static List<InternalFragmentIon> Run(
            string psmTsvPath,
            string rawFileFolder,
            string outputDirectory,
            double defaultCollisionEnergy = double.NaN)
        {
            // Validate inputs
            if (!File.Exists(psmTsvPath))
                throw new FileNotFoundException($"PSM TSV file not found: {psmTsvPath}");
            if (!Directory.Exists(rawFileFolder))
                throw new DirectoryNotFoundException($"Raw file folder not found: {rawFileFolder}");
            if (!Directory.Exists(outputDirectory))
                Directory.CreateDirectory(outputDirectory);

            // Load PSMs
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);

            // Load raw files
            var msDataFiles = LoadRawFiles(rawFileFolder, psms);

            // Extract internal fragments
            var allInternalIons = new List<InternalFragmentIon>();
            var psmsByFile = psms.GroupBy(p => p.FileNameWithoutExtension).ToList();

            foreach (var fileGroup in psmsByFile)
            {
                string fileName = fileGroup.Key;
                var filePsms = fileGroup.ToList();

                if (!msDataFiles.TryGetValue(fileName, out var msDataFile))
                    continue;

                var internalIons = InternalFragmentFeatureExtractor.ExtractFromPsms(
                    filePsms, msDataFile, defaultCollisionEnergy);
                allInternalIons.AddRange(internalIons);
            }

            // Write output
            string outputPath = Path.Combine(outputDirectory, "InternalFragmentIons.tsv");
            InternalFragmentTsvWriter.WriteToTsv(allInternalIons, outputPath);

            return allInternalIons;
        }

        private static Dictionary<string, MsDataFile> LoadRawFiles(
            string rawFileFolder,
            List<PsmFromTsv> psms)
        {
            var supportedExtensions = new[] { ".raw", ".mzml", ".mgf" };
            var rawFiles = supportedExtensions
                .SelectMany(ext => Directory.GetFiles(rawFileFolder, $"*{ext}", SearchOption.TopDirectoryOnly))
                .ToList();

            var neededFiles = psms
                .Select(p => p.FileNameWithoutExtension)
                .Distinct()
                .ToHashSet(StringComparer.OrdinalIgnoreCase);

            var msDataFiles = new Dictionary<string, MsDataFile>(StringComparer.OrdinalIgnoreCase);

            foreach (var rawFile in rawFiles)
            {
                string fileNameWithoutExt = Path.GetFileNameWithoutExtension(rawFile);

                if (!neededFiles.Contains(fileNameWithoutExt))
                    continue;

                try
                {
                    string extension = Path.GetExtension(rawFile).ToLowerInvariant();
                    MsDataFile msDataFile = extension switch
                    {
                        ".mzml" => Mzml.LoadAllStaticData(rawFile),
                        ".raw" => ThermoRawFileReader.LoadAllStaticData(rawFile),
                        ".mgf" => Mgf.LoadAllStaticData(rawFile),
                        _ => throw new NotSupportedException($"Unsupported file format: {extension}")
                    };

                    msDataFiles[fileNameWithoutExt] = msDataFile;
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"Warning: Could not load {rawFile}: {ex.Message}");
                }
            }

            return msDataFiles;
        }
    }
}