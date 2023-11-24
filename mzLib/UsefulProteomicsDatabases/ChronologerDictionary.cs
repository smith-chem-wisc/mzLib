using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class ChronologerDictionary
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
            };


            if (dict == TypeOfDictionary.Unimod)
            {
                int aaCount = 20;

                var mods =
                    UnimodLoader.ReadMods(@"F:\Research\Data\unimod.xml").ToList();

                var groupedModsByOriginalID =
                    mods.GroupBy(x => x.Target.ToString()).ToList();

                foreach (var target in groupedModsByOriginalID)
                {
                    foreach (var mod in target)
                    {
                        aaCount = aaCount + 1;
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
                    //_residueWithModToTensorInt.Add(('C',null),23);//'S - carbamidomethylcysteine
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
