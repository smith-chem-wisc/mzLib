using Proteomics.Fragmentation;
using System;

namespace Proteomics.ProteolyticDigestion
{
    public class DigestionParams
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this("trypsin")
        {
        }

        public DigestionParams(string protease = "trypsin", int maxMissedCleavages = 2, int minPeptideLength = 7, int maxPeptideLength = int.MaxValue,
            int maxModificationIsoforms = 1024, InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
            int maxModsForPeptides = 2, bool semiSpecificDigestion = false, bool nonSpecificDigestion = false)
        {
            Protease = ProteaseDictionary.Dictionary[protease];
            MaxMissedCleavages = maxMissedCleavages;
            MinPeptideLength = minPeptideLength;
            MaxPeptideLength = maxPeptideLength;
            MaxModificationIsoforms = maxModificationIsoforms;
            InitiatorMethionineBehavior = initiatorMethionineBehavior;
            MaxModsForPeptide = maxModsForPeptides;
            SemiSpecificDigestion = semiSpecificDigestion;
            NonSpecificDigestion = nonSpecificDigestion;
        }
        
        public int MaxMissedCleavages { get; private set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MinPeptideLength { get; private set; }
        public int MaxPeptideLength { get; private set; }
        public int MaxModificationIsoforms { get; private set; }
        public int MaxModsForPeptide { get; private set; }
        public Protease Protease { get; private set; }
        public bool SemiSpecificDigestion { get; private set; } //for semispecific searching of proteases
        public bool NonSpecificDigestion { get; private set; } //for nonspecific searching of proteases (where the digestion is non-specific, but the protease is used for FDR)

        public override bool Equals(object obj)
        {
            DigestionParams a = obj as DigestionParams;

            return a != null
                && MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && MinPeptideLength.Equals(a.MinPeptideLength)
                && MaxPeptideLength.Equals(a.MaxPeptideLength)
                && InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && Protease.Equals(a.Protease)
                && SemiSpecificDigestion.Equals(a.SemiSpecificDigestion)
                && NonSpecificDigestion.Equals(a.NonSpecificDigestion);
        }
        
        public override int GetHashCode()
        {
            return
                MaxMissedCleavages.GetHashCode()
                ^ InitiatorMethionineBehavior.GetHashCode()
                ^ MaxModificationIsoforms.GetHashCode()
                ^ MaxModsForPeptide.GetHashCode();
        }

        public override string ToString()
        {
            return MaxMissedCleavages + "," + InitiatorMethionineBehavior + "," + MinPeptideLength + "," + MaxPeptideLength + "," 
                + MaxModificationIsoforms + "," + MaxModsForPeptide + "," + Protease.Name + "," + SemiSpecificDigestion + ","
                + NonSpecificDigestion;
        }

        /// <summary>
        /// Creates a DigestionParams object from string. Used after deserializing a PeptideWithSetModifications
        /// </summary>
        public static DigestionParams FromString(string str)
        {
            string[] split = str.Split(',');
            return new DigestionParams(
                protease: split[6],
                maxMissedCleavages: int.Parse(split[0]),
                minPeptideLength: int.Parse(split[2]),
                maxPeptideLength: int.Parse(split[3]),
                maxModificationIsoforms: int.Parse(split[4]),
                initiatorMethionineBehavior: (InitiatorMethionineBehavior)Enum.Parse(typeof(InitiatorMethionineBehavior), split[1]),
                maxModsForPeptides: int.Parse(split[5]),
                semiSpecificDigestion: bool.Parse(split[7]),
                nonSpecificDigestion: bool.Parse(split[8]));
        }
    }
}