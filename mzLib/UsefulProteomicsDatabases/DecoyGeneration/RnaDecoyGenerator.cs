using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Omics.Modifications;
using Transcriptomics;
using Omics.BioPolymer;

namespace UsefulProteomicsDatabases
{
    public static class RnaDecoyGenerator
    {
        public static List<T> GenerateDecoys<T>(List<T> nucleicAcids, DecoyType decoyType, int maxThreads = -1, string decoyIdentifier = "DECOY") where T : INucleicAcid
        {
            switch (decoyType)
            {
                case DecoyType.None:
                    return new List<T>();
                case DecoyType.Reverse:
                    return GenerateReverseDecoys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Slide:
                    return GenerateSlidedDecoys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Shuffle:
                    return GenerateShuffledDeocys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Random:
                default:
                    throw new ArgumentOutOfRangeException(nameof(decoyType), decoyType, null);
            }
        }

        /// <summary>
        /// Reverse decoys: sequence reversed, 3' terminus retained chemically (termini objects preserved),
        /// modifications & variant-specific modifications follow their original nucleotide.
        /// Each modification is cloned with a motif matching the nucleotide at its new (reversed) coordinate
        /// to avoid motif/base mismatch filtering during RNA construction.
        /// </summary>
        private static List<T> GenerateReverseDecoys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            List<T> decoyNucleicAcids = new();
            Parallel.ForEach(nucleicAcids, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, nucleicAcid =>
            {
                string originalSeq = nucleicAcid.BaseSequence;
                int L = originalSeq.Length;

                // Reverse sequence characters
                string reverseSequence = new string(originalSeq.Reverse().ToArray());

                // Map original 1-based index -> reversed 1-based index
                Dictionary<int, int> indexMapping = new(L);
                for (int i = 1; i <= L; i++)
                {
                    indexMapping[i] = L - i + 1;
                }

                // Helper: try to clone a modification for a specific nucleotide.
                // If cloning fails (constructor signature mismatches), we fall back to the original modification instance.
                static Modification CloneForBase(Modification mod, char nucleotide)
                {
                    if (!ModificationMotif.TryGetMotif(nucleotide.ToString(), out var motif))
                    {
                        // Fallback: reuse existing motif (may be null)
                        motif = mod.Target;
                    }

                    try
                    {
                        // Prefer the most common simple constructor.
                        // Many test-created modifications use a short signature:
                        // (originalId, something?, modificationType, something?, motif, locationRestriction, formula?)
                        // We only preserve OriginalId, ModificationType (if available), motif, and location restriction when accessible.
                        string originalId = mod.OriginalId ?? mod.IdWithMotif ?? "";
                        string modificationType = mod.ModificationType ?? "Cloned";
                        string locationRestriction = mod.LocationRestriction ?? "Anywhere.";

                        // Attempt to keep formula & masses if available
                        var formula = mod.ChemicalFormula; // may be null
                        if (formula != null)
                        {
                            return new Modification(
                                _originalId: originalId,
                                _modificationType: modificationType,
                                _target: motif,
                                _locationRestriction: locationRestriction,
                                _chemicalFormula: formula);
                        }
                        // Fallback minimal
                        return new Modification(
                            _originalId: originalId,
                            _modificationType: modificationType,
                            _target: motif,
                            _locationRestriction: locationRestriction);
                    }
                    catch
                    {
                        // Fallback: return original if construction path unknown
                        return mod;
                    }
                }

                // Reverse base-level modifications by cloning for the nucleotide that moves.
                var reverseModifications = new Dictionary<int, List<Modification>>();
                foreach (var kvp in nucleicAcid.OneBasedPossibleLocalizedModifications)
                {
                    int originalIndex = kvp.Key;
                    int reversedIndex = indexMapping[originalIndex];
                    char nucleotide = originalSeq[originalIndex - 1];

                    var clonedList = new List<Modification>(kvp.Value.Count);
                    foreach (var m in kvp.Value)
                    {
                        clonedList.Add(CloneForBase(m, nucleotide));
                    }
                    reverseModifications[reversedIndex] = clonedList;
                }

                List<TruncationProduct> reverseTruncs = new();
                List<SequenceVariation> reverseVariations = new();
                List<SequenceVariation> reverseAppliedVariations = new();

                if (nucleicAcid is IHasSequenceVariants variantContaining)
                {
                    static void Normalize(ref int a, ref int b)
                    {
                        if (a > b) (a, b) = (b, a);
                    }

                    SequenceVariation ReverseVariant(SequenceVariation v)
                    {
                        int rb = indexMapping[v.OneBasedBeginPosition];
                        int re = indexMapping[v.OneBasedEndPosition];
                        Normalize(ref rb, ref re);

                        // Reverse variant-specific modifications
                        Dictionary<int, List<Modification>> reversedVariantMods = null;
                        if (v.OneBasedModifications != null && v.OneBasedModifications.Count > 0)
                        {
                            reversedVariantMods = new Dictionary<int, List<Modification>>(v.OneBasedModifications.Count);
                            foreach (var modKvp in v.OneBasedModifications)
                            {
                                int revKey = indexMapping[modKvp.Key];
                                char baseChar = originalSeq[modKvp.Key - 1];
                                var cloned = modKvp.Value.Select(m => CloneForBase(m, baseChar)).ToList();
                                reversedVariantMods[revKey] = cloned;
                            }
                        }

                        return new SequenceVariation(
                            rb,
                            re,
                            v.OriginalSequence,
                            v.VariantSequence,
                            v.Description,
                            v.VariantCallFormatData?.Description,
                            reversedVariantMods);
                    }

                    foreach (var v in variantContaining.AppliedSequenceVariations)
                    {
                        reverseAppliedVariations.Add(ReverseVariant(v));
                    }

                    foreach (var v in variantContaining.SequenceVariations)
                    {
                        reverseVariations.Add(ReverseVariant(v));
                    }

                    // Reverse truncations
                    foreach (var t in variantContaining.TruncationProducts)
                    {
                        if (t.OneBasedBeginPosition.HasValue && t.OneBasedEndPosition.HasValue)
                        {
                            int rb = indexMapping[t.OneBasedEndPosition.Value];
                            int re = indexMapping[t.OneBasedBeginPosition.Value];
                            Normalize(ref rb, ref re);
                            reverseTruncs.Add(new TruncationProduct(rb, re, $"{decoyIdentifier} {t.Type}"));
                        }
                    }
                }

                // Construct decoy
                T newNucleicAcid = nucleicAcid.CreateNew(
                    reverseSequence,
                    reverseModifications,
                    isDecoy: true,
                    truncationProducts: reverseTruncs,
                    sequenceVariations: reverseVariations,
                    appliedSequenceVariations: reverseAppliedVariations,
                    decoyIdentifier: decoyIdentifier);

                lock (decoyNucleicAcids)
                {
                    decoyNucleicAcids.Add(newNucleicAcid);
                }
            });
            return decoyNucleicAcids;
        }

        private static List<T> GenerateSlidedDecoys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }

        private static List<T> GenerateShuffledDeocys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }
    }
}
