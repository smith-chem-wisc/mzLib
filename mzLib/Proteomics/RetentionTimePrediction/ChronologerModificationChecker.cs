using System.Collections.Generic;
using System.Linq;
using System.Text;
using Chromatography.RetentionTimePrediction.Util;

namespace Proteomics.RetentionTimePrediction
{
    /// <summary>
    /// Checks modification compatibility for Chronologer predictor.
    /// Lives in Proteomics layer since it depends on Modification class.
    /// </summary>
    public class ChronologerModificationChecker : IModificationCompatibilityChecker
    {
        private static readonly HashSet<string> SupportedModifications = new()
        {
            "Carbamidomethyl on C",
            "Oxidation on M",
            "Glu to PyroGlu",
            "Gln to PyroGlu",
            "Phosphorylation on S",
            "Phosphorylation on T",
            "Phosphorylation on Y",
            "Acetylation on K",
            "Succinylation on K",
            "Ubiquitination on K",
            "Methylation on K",
            "Dimethylation on K",
            "Trimethylation on K",
            "Methylation on R",
            "Dimethylation on R",
            "GlyGly on K"
        };

        public bool AreModificationsCompatible(string fullSequence, out List<string>? incompatibleModIds)
        {
            // Parse modifications from full sequence
            var mods = ExtractModificationIds(fullSequence);
            
            incompatibleModIds = mods
                .Where(modId => !SupportedModifications.Contains(modId))
                .ToList();

            return incompatibleModIds.Count == 0;
        }

        public string FilterIncompatibleModifications(string fullSequence)
        {
            // Remove incompatible modification annotations from full sequence
            var sb = new StringBuilder();
            int i = 0;
            
            while (i < fullSequence.Length)
            {
                if (fullSequence[i] == '[')
                {
                    // Found modification
                    int closeIdx = fullSequence.IndexOf(']', i);
                    if (closeIdx != -1)
                    {
                        string modAnnotation = fullSequence.Substring(i + 1, closeIdx - i - 1);
                        string modId = ExtractModIdFromAnnotation(modAnnotation);
                        
                        if (SupportedModifications.Contains(modId))
                        {
                            // Keep this modification
                            sb.Append(fullSequence, i, closeIdx - i + 1);
                        }
                        // else skip (filter out)
                        
                        i = closeIdx + 1;
                        continue;
                    }
                }
                
                sb.Append(fullSequence[i]);
                i++;
            }
            
            return sb.ToString();
        }

        public string FilterIncompatibleModificationsFromMassShifts(
            string originalMassShiftSequence,
            string fullSequence)
        {
            // Extract incompatible mod positions from full sequence
            var incompatiblePositions = GetIncompatibleModificationPositions(fullSequence);
            
            // Remove mass shifts at those positions
            return RemoveMassShiftsAtPositions(originalMassShiftSequence, incompatiblePositions);
        }

        private List<string> ExtractModificationIds(string fullSequence)
        {
            var mods = new List<string>();
            int i = 0;
            
            while (i < fullSequence.Length)
            {
                if (fullSequence[i] == '[')
                {
                    int closeIdx = fullSequence.IndexOf(']', i);
                    if (closeIdx != -1)
                    {
                        string modAnnotation = fullSequence.Substring(i + 1, closeIdx - i - 1);
                        string modId = ExtractModIdFromAnnotation(modAnnotation);
                        mods.Add(modId);
                        i = closeIdx + 1;
                        continue;
                    }
                }
                i++;
            }
            
            return mods;
        }

        private string ExtractModIdFromAnnotation(string annotation)
        {
            // Format: "Variable:Oxidation on M" or "Fixed:Carbamidomethyl on C"
            int colonIdx = annotation.IndexOf(':');
            if (colonIdx >= 0)
            {
                return annotation.Substring(colonIdx + 1);
            }
            return annotation;
        }

        private HashSet<int> GetIncompatibleModificationPositions(string fullSequence)
        {
            var positions = new HashSet<int>();
            int residuePosition = 0;
            int i = 0;
            
            while (i < fullSequence.Length)
            {
                if (fullSequence[i] == '[')
                {
                    int closeIdx = fullSequence.IndexOf(']', i);
                    if (closeIdx != -1)
                    {
                        string modAnnotation = fullSequence.Substring(i + 1, closeIdx - i - 1);
                        string modId = ExtractModIdFromAnnotation(modAnnotation);
                        
                        if (!SupportedModifications.Contains(modId))
                        {
                            positions.Add(residuePosition);
                        }
                        
                        i = closeIdx + 1;
                        continue;
                    }
                }
                
                if (char.IsLetter(fullSequence[i]))
                {
                    residuePosition++;
                }
                
                i++;
            }
            
            return positions;
        }

        private string RemoveMassShiftsAtPositions(string massShiftSeq, HashSet<int> positions)
        {
            var sb = new StringBuilder();
            int residuePosition = 0;
            int i = 0;
            
            while (i < massShiftSeq.Length)
            {
                if (massShiftSeq[i] == '[')
                {
                    int closeIdx = massShiftSeq.IndexOf(']', i);
                    if (closeIdx != -1)
                    {
                        if (!positions.Contains(residuePosition))
                        {
                            // Keep this mass shift
                            sb.Append(massShiftSeq, i, closeIdx - i + 1);
                        }
                        // else skip (filter out)
                        
                        i = closeIdx + 1;
                        continue;
                    }
                }
                
                sb.Append(massShiftSeq[i]);
                
                if (char.IsLetter(massShiftSeq[i]))
                {
                    residuePosition++;
                }
                
                i++;
            }
            
            return sb.ToString();
        }
    }
}