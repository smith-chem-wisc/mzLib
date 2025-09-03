using System.Text;
using System.Text.RegularExpressions;
using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;

namespace Omics
{
    /// <summary>
    /// Interface for modified and unmodified precursor ions
    /// </summary>
    /// <remarks>
    /// Proteins -> PeptideWithSetModifications : ProteolyticPeptide
    /// Nucleic Acids -> OligoWithSetMods : NucleolyticOligo
    /// </remarks>
    public interface IBioPolymerWithSetMods : IHasChemicalFormula, IEquatable<IBioPolymerWithSetMods>
    {
        string BaseSequence { get; }
        string FullSequence { get; }
        double MostAbundantMonoisotopicMass { get; }
        string SequenceWithChemicalFormulas { get; }
        int OneBasedStartResidue { get; }
        int OneBasedEndResidue { get; }
        int MissedCleavages { get; }
        /// <summary>
        /// Description of where the BioPolymerWithSetMods originated from examples include
        /// Top-down truncation: full-length proteoform C-terminal digestion truncation
        /// Top-down truncation: DECOY full-length proteoform N-terminal digestion truncation
        /// Bottom-up search: full
        /// Bottom-up search: DECOY full
        /// Bottom-up search : chain(49-597) start
        /// </summary>
        string Description { get; }
        CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; }
        char PreviousResidue { get; }
        char NextResidue { get; }
        IDigestionParams DigestionParams { get; }
        Dictionary<int, Modification> AllModsOneIsNterminus { get; }
        int NumMods { get; }
        int NumFixedMods { get; }
        int NumVariableMods { get; }
        int Length { get; }
        char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];
        IBioPolymer Parent { get; }

        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
            List<Product> products);

        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<Product> products);

        /// <summary>
        /// Outputs a duplicate IBioPolymerWithSetMods with a localized mass shift, replacing a modification when present
        /// <remarks>
        /// Used to localize an unknown mass shift in the MetaMorpheus Localization Engine
        /// </remarks>
        /// </summary>
        /// <param name="indexOfMass">The index of the modification in the AllModOneIsNTerminus Dictionary - 2 (idk why -2)</param>
        /// <param name="massToLocalize">The mass to add to the BioPolymer</param>
        /// <returns></returns>
        public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize);

        public static string GetBaseSequenceFromFullSequence(string fullSequence, char modStartDelimiter = '[', char modEndDelimiter = ']', char cTerminusDelimiter = '-')
        {
            StringBuilder sb = new StringBuilder();
            int bracketCount = 0;
            foreach (char c in fullSequence)
            {
                if (c == modStartDelimiter)
                {
                    bracketCount++;
                }
                else if (c == modEndDelimiter)
                {
                    bracketCount--;
                }
                else if (c == cTerminusDelimiter)
                {
                    continue;
                }
                else if (bracketCount == 0)
                {
                    sb.Append(c);
                }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Returns a list of modifications and their OneBased index from a full sequence
        /// </summary>
        /// <param name="fullSequence">Full sequence</param>
        /// <param name="allModsKnown">All known modifications</param>
        /// <returns></returns>
        /// <exception cref="MzLibUtil.MzLibException">When a full sequence is not in the correct format or a mod is not found in the allModsKnown dictionary</exception>
        public static Dictionary<int, Modification> GetModificationDictionaryFromFullSequence(string fullSequence,
            Dictionary<string, Modification> allModsKnown)
        {
            var allModsOneIsNterminus = new Dictionary<int, Modification>();
            var baseSequence = GetBaseSequenceFromFullSequence(fullSequence);
            int currentModStart = 0;
            int currentModificationLocation = 1;
            bool currentlyReadingMod = false;
            int bracketCount = 0;

            for (int r = 0; r < fullSequence.Length; r++)
            {
                char c = fullSequence[r];
                if (c == '[')
                {
                    currentlyReadingMod = true;
                    if (bracketCount == 0)
                    {
                        currentModStart = r + 1;
                    }
                    bracketCount++;
                }
                else if (c == ']')
                {
                    string modId = null;
                    bracketCount--;
                    if (bracketCount == 0)
                    {
                        try
                        {
                            //remove the beginning section (e.g. "Fixed", "Variable", "Uniprot")
                            string modString = fullSequence.Substring(currentModStart, r - currentModStart);
                            int splitIndex = modString.IndexOf(':');
                            string modType = modString.Substring(0, splitIndex);
                            modId = modString.Substring(splitIndex + 1, modString.Length - splitIndex - 1);
                        }
                        catch (Exception e)
                        {
                            throw new MzLibUtil.MzLibException(
                                "Error while trying to parse string into peptide: " + e.Message, e);

                        }
                        if (!allModsKnown.TryGetValue(modId, out var mod))
                        {
                            throw new MzLibUtil.MzLibException(
                                "Could not find modification while reading string: " + fullSequence);
                        }
                        if (mod.LocationRestriction.Contains("C-terminal.") && r == fullSequence.Length - 1)
                        {
                            currentModificationLocation = baseSequence.Length + 2;
                        }
                        allModsOneIsNterminus.Add(currentModificationLocation, mod);
                        currentlyReadingMod = false;
                    }
                }
                else if (!currentlyReadingMod)
                {
                    currentModificationLocation++;
                }
                //else do nothing
            }

            return allModsOneIsNterminus;
        }

        /// <summary>
        /// Returns a list of modifications from a full sequence
        /// </summary>
        /// <param name="fullSequence">Full sequence</param>
        /// <param name="allModsKnown">All known modifications</param>
        /// <returns></returns>
        public static List<Modification> GetModificationsFromFullSequence(string fullSequence,
            Dictionary<string, Modification> allModsKnown) => [.. GetModificationDictionaryFromFullSequence(fullSequence, allModsKnown).Values];

        /// <summary>
        /// Determines the full sequence of a BioPolymerWithSetMods from its base sequence and modifications
        /// </summary>
        public static string DetermineFullSequence(string baseSequence, IDictionary<int, Modification> allModsOneIsNterminus)
        {
            // start string builder with initial capacity to avoid resizing costs. 
            var subSequence = new StringBuilder(baseSequence.Length + allModsOneIsNterminus.Count * 60);

            // modification on peptide N-terminus
            if (allModsOneIsNterminus.TryGetValue(1, out Modification? mod))
            {
                subSequence.Append($"[{mod.ModificationType}:{mod.IdWithMotif}]");
            }

            for (int r = 0; r < baseSequence.Length; r++)
            {
                subSequence.Append(baseSequence[r]);

                // modification on this residue
                if (allModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subSequence.Append($"[{mod.ModificationType}:{mod.IdWithMotif}]");
                }
            }

            // modification on peptide C-terminus
            if (allModsOneIsNterminus.TryGetValue(baseSequence.Length + 2, out mod))
            {
                subSequence.Append($"-[{mod.ModificationType}:{mod.IdWithMotif}]");
            }

            return subSequence.ToString();
        }

        public static string ParseSubstitutedFullSequence(string fullSequenceWithSubstitution)
        {
            string pattern = @"\[\d+\+?\s*nucleotide substitution:\s*([A-Z])->([A-Z]) on ([A-Z])\]";
            string parsedSequence = fullSequenceWithSubstitution;
            var match = Regex.Match(parsedSequence, pattern);
            while (match.Success)
            {
                string original = match.Groups[1].Value; // Z (original)
                string sub = match.Groups[2].Value;   // Y (substitute)
                int patternIndex = match.Index;
                int replaceIndex = parsedSequence.LastIndexOf(original, patternIndex);
                if (replaceIndex != -1)
                {
                    // Replace the first occurrence of Z before the pattern with Y
                    parsedSequence = parsedSequence.Remove(replaceIndex, 1).Insert(replaceIndex, sub);
                }
                // Remove the pattern
                parsedSequence = parsedSequence.Remove(patternIndex, match.Length);
                match = Regex.Match(parsedSequence, pattern); // Find the next match
            }
            return parsedSequence;
        }
    }
}
