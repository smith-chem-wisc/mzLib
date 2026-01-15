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
            var assembly = Assembly.GetExecutingAssembly();

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
                    string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseProteaseLines(lines.Skip(1).ToArray(), proteaseMods);
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
            myLines = myLines.Skip(1).ToArray();
            return ParseProteaseLines(myLines, proteaseMods);
        }

        /// <summary>
        /// Parses protease definitions from TSV-formatted lines.
        /// </summary>
        /// <param name="lines">Lines from the proteases file (header already skipped).</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        private static Dictionary<string, Protease> ParseProteaseLines(string[] lines, List<Modification> proteaseMods)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();

            foreach (string line in lines)
            {
                if (string.IsNullOrWhiteSpace(line))
                {
                    continue; // skip empty lines
                }

                string[] fields = line.Split('\t');
                string name = fields[0];
                List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(fields[1]);
                var cleavageSpecificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true);
                string psiMsAccessionNumber = fields[5];
                string psiMsName = fields[6];
                string proteaseModDetails = fields[8]; // name of the modification associated with proteolytic cleavage

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

                // Modification was specified but not found - throw after adding the protease without the mod
                throw new MzLibException($"{proteaseModDetails} is not a valid modification");
            }

            // No modification required
            return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
        }
    }
}