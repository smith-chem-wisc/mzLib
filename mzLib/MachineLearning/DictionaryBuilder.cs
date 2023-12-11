using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// This class provides static methods to build dictionaries for different models.
    /// </summary>
    public static class DictionaryBuilder
    {
        public static Dictionary<(char, string), int>
            GetChronologerDictionary(TypeOfDictionary dict)
        {
            var dictionary = new Dictionary<(char, string), int>()
            {
                { ('A', ""), 1  }, //'Alanine
                { ('C', ""), 2  }, //'Cysteine
                { ('D', ""), 3  }, //'Aspartate
                { ('E', ""), 4  }, //'Glutamate
                { ('F', ""), 5  }, //'Phenylalaline
                { ('G', ""), 6  }, //'Glycine
                { ('H', ""), 7  }, //'Histidine
                { ('I', ""), 8  }, //'Isoleucine
                { ('K', ""), 9  }, //'Lysine
                { ('L', ""), 10 }, //'Leucine
                { ('M', ""), 11 }, //'Methionine
                { ('N', ""), 12 }, //'Asparagine
                { ('P', ""), 13 }, //'Proline
                { ('Q', ""), 14 }, //'Glutamine
                { ('R', ""), 15 }, //'Argenine
                { ('S', ""), 16 }, //'Serine
                { ('T', ""), 17 }, //'Threonine
                { ('V', ""), 18 }, //'Valine
                { ('W', ""), 19 }, //'Tryptophane
                { ('Y', ""), 20 }, //'Tyrosine
                {('C', "C-Terminus"), 38},
                {('N', "N-Terminus"), 44}
            };

            //todo: change numbers to accession numbers from DB
            if (dict == TypeOfDictionary.Unimod)
            {
                int aaCount = 20;

                var mods =
                    Loaders.LoadUnimod(@"F:\Research\Data\unimod.xml").ToList();

                var groupedModsByOriginalID =
                    mods.GroupBy(x => x.Target.ToString()).ToList();

                foreach (var target in groupedModsByOriginalID)
                {
                    foreach (var mod in target)
                    {
                        aaCount = aaCount + 1;

                        //C-Terminus and N-Terminus
                        if(aaCount == 38 || aaCount == 44)
                        {
                            aaCount = aaCount + 1;
                        }

                        if (!dictionary.ContainsKey((Char.Parse(target.Key), mod.IdWithMotif)))
                            dictionary.Add((Char.Parse(target.Key), mod.IdWithMotif), aaCount);
                    }
                }

                return dictionary;
            }

            if (dict == TypeOfDictionary.Chronologer)
            {
                var chronologerModsDict = new Dictionary<(char, string), int>()
                {
                    { ('C', "Carbamidomethyl on C"), 21 }, //'Carbamidomethyl
                    { ('M', "Oxidation on M"), 22 }, //'Oxidized
                    { ('C', "Pyro-carbamidomethyl on C"), 23 },//'S - carbamidomethylcysteine
                    { ('E', "Glu to PyroGlu"), 24 }, //'Pyroglutamate
                    { ('S', "Phosphorylation on S"), 25 }, //'Phosphoserine
                    { ('T', "Phosphorylation on T"), 26 }, //'Phosphothreonine
                    { ('Y', "Phosphorylation on Y"), 27 }, //'Phosphotyrosine
                    { ('K', "Accetylation on K"), 28 }, //'Acetylated
                    { ('K', "Succinylation on K"), 29 }, //'Succinylated
                    { ('K', "Ubiquitination on K"), 30 }, //'Ubiquitinated
                    { ('K', "Methylation on K"), 31 }, //'Monomethyl
                    { ('K', "Dimethylation on K"), 32 }, //'Dimethyl
                    { ('K', "Trimethylation on K"), 33 }, //'Trimethyl
                    { ('R', "Methylation on R"), 34 }, //'Monomethyl
                    { ('R', "Dimethylation on R"), 35 }, //'Dimethyl
                    { ('K', "TMT on K"), 36 }, //TMT0-modified lysisne
                    { ('K', "TMT6plex on K"), 37 }, //tmt10-modified lysine

                };

                foreach (var item in chronologerModsDict)
                {
                    if (!dictionary.ContainsKey(item.Key))
                        dictionary.Add(item.Key, item.Value);
                }
            }
            
            if (dict == TypeOfDictionary.CanonicalAA)
            {
                return dictionary;
            }

            return dictionary;
        }

        public enum TypeOfDictionary
        {
            CanonicalAA,
            Chronologer,
            Unimod
        }
    }
}
