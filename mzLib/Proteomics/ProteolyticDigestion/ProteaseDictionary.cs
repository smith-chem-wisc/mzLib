using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;

namespace Proteomics.ProteolyticDigestion
{
    public static class ProteaseDictionary
    {
        private const string EmbeddedResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";

        static ProteaseDictionary()
        {
            // Load from embedded resource (no protease modifications in static initialization)
            Dictionary = LoadProteaseDictionary(proteaseMods: null);
        }

        public static Dictionary<string, Protease> Dictionary { get; set; }

        /// <summary>
        /// Loads the default proteases from the embedded resource.
        /// </summary>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        public static Dictionary<string, Protease> LoadProteaseDictionary(List<Modification> proteaseMods)
        {
            var assembly = typeof(ProteaseDictionary).Assembly;

            using (var stream = assembly.GetManifestResourceStream(EmbeddedResourceName))
            {
                if (stream == null)
                {
                    throw new MzLibException(
                        $"Could not find embedded resource '{EmbeddedResourceName}'. " +
                        $"Available resources: {string.Join(", ", assembly.GetManifestResourceNames())}");
                }

                using (var reader = new StreamReader(stream))
                {
                    string fileContent = reader.ReadToEnd();
                    // RemoveEmptyEntries skips blank lines and lines with only whitespace,
                    // which is the desired behavior for TSV parsing
                    string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseProteaseLines(lines, proteaseMods);
                }
            }
        }

        /// <summary>
        /// Loads proteases from an external file path. Useful for loading custom user-defined proteases.
        /// </summary>
        /// <param name="path">Path to the proteases.tsv file.</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        public static Dictionary<string, Protease> LoadProteaseDictionary(string path, List<Modification> proteaseMods = null)
        {
            string[] myLines = File.ReadAllLines(path);
            return ParseProteaseLines(myLines, proteaseMods);
        }

        /// <summary>
        /// Parses protease definitions from TSV-formatted lines.
        /// Lines starting with '#' are treated as comments and skipped.
        /// The header line (starting with "Name") is also skipped.
        /// </summary>
        /// <param name="lines">Lines from the proteases file.</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        private static Dictionary<string, Protease> ParseProteaseLines(string[] lines, List<Modification> proteaseMods)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();

            foreach (string line in lines)
            {
                // Skip empty lines, comment lines, and the header line
                // Trim to handle potential BOM or leading/trailing whitespace
                string trimmedLine = line.Trim().TrimStart('\uFEFF'); // \uFEFF is the UTF-8 BOM character
                
                if (string.IsNullOrWhiteSpace(trimmedLine) || 
                    trimmedLine.StartsWith("#") || 
                    trimmedLine.StartsWith("Name\t", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }

                string[] fields = trimmedLine.Split('\t');
                
                if (fields.Length < 3)
                {
                    throw new MzLibException($"Protease definition line has insufficient fields (expected at least 3, got {fields.Length}): {line}");
                }

                // Expected TSV columns (see proteases.tsv header for full documentation):
                // [0] Name                  - Required: Unique protease identifier
                // [1] Motif                 - Required: Cleavage motif syntax (e.g., "K[P]|,R[P]|")
                // [2] Specificity           - Required: full, semi, none, SingleN, or SingleC
                // [3] PSI-MS Accession      - Optional: Standard identifier (e.g., MS:1001313)
                // [4] PSI-MS Name           - Optional: Standard name from PSI-MS ontology
                // [5] Cleavage Modification - Optional: Modification applied at cleavage site
                
                // Trim each field to handle whitespace and ensure empty fields are truly empty
                string name = fields[0].Trim();
                string motifField = fields[1].Trim();
                string specificityField = fields[2].Trim();
                string psiMsAccessionNumber = fields.Length > 3 ? fields[3].Trim() : string.Empty;
                string psiMsName = fields.Length > 4 ? fields[4].Trim() : string.Empty;
                string proteaseModDetails = fields.Length > 5 ? fields[5].Trim() : string.Empty;

                List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(motifField);
                var cleavageSpecificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), specificityField, true);

                Protease protease = CreateProtease(
                    name, motifList, cleavageSpecificity, psiMsAccessionNumber, psiMsName,
                    proteaseModDetails, proteaseMods);

                if (dict.ContainsKey(protease.Name))
                {
                    throw new MzLibException($"More than one protease named {protease.Name} exists");
                }

                dict.Add(protease.Name, protease);
            }

            return dict;
        }

        /// <summary>
        /// Creates a Protease object, optionally with an associated modification.
        /// </summary>
        private static Protease CreateProtease(
            string name,
            List<DigestionMotif> motifList,
            CleavageSpecificity cleavageSpecificity,
            string psiMsAccessionNumber,
            string psiMsName,
            string proteaseModDetails,
            List<Modification> proteaseMods)
        {
            // If this protease has an associated modification, look it up
            if (!string.IsNullOrEmpty(proteaseModDetails) && proteaseMods != null)
            {
                Modification proteaseModification = proteaseMods
                    .FirstOrDefault(p => p.IdWithMotif == proteaseModDetails);

                if (proteaseModification != null)
                {
                    return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList, proteaseModification);
                }

                // Modification was specified but not found in the provided list
                throw new MzLibException($"{proteaseModDetails} is not a valid modification");
            }

            // No modification required
            return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
        }
    }
}