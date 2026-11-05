using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Random proteins, modifications and digestion parameters for the equivalence tests of digestion: proteins carry
    /// localized modifications (some of them terminal, some sharing a residue), sequences occasionally hold residues of
    /// unknown mass, and the parameters cover full, semi and top-down digestion, every initiator-methionine behaviour,
    /// and isoform caps low enough to cut enumeration short.
    /// </summary>
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal static class RandomDigestionInputs
    {
        // K and R twice, so tryptic products are of ordinary length; X and B rarely, for residues of unknown mass
        private const string Residues = "ACDEFGHIKLMNPQRSTVWYKRACDEFGHIKLMNPQRSTVWYKRUXB";

        internal static readonly Modification CarbamidomethylC = New("Carbamidomethyl", "C", "Anywhere.", 57.021464);
        internal static readonly Modification OxidationM = New("Oxidation", "M", "Anywhere.", 15.994915);
        internal static readonly Modification PhosphoS = New("Phospho", "S", "Anywhere.", 79.966331);
        internal static readonly Modification AcetylNTerm = New("Acetyl", "X", "N-terminal.", 42.010565);

        internal static readonly List<Modification> Localized = new()
        {
            PhosphoS,
            New("Phospho", "T", "Anywhere.", 79.966331),
            New("Phospho", "Y", "Anywhere.", 79.966331),
            New("Acetyl", "K", "Anywhere.", 42.010565),
            New("Methyl", "K", "Anywhere.", 14.01565),
            New("Methyl", "R", "Anywhere.", 14.01565),
            OxidationM,
            New("Deamidation", "N", "Anywhere.", 0.984016),
            New("Hydroxylation", "P", "Anywhere.", 15.994915),
            AcetylNTerm,
            New("Amidation", "X", "C-terminal.", -0.984016),
            New("Formyl", "X", "Peptide N-terminal.", 27.994915),
            New("Methylester", "X", "Peptide C-terminal.", 14.01565),
        };

        internal static Modification New(string id, string motif, string location, double mass)
        {
            ModificationMotif.TryGetMotif(motif, out ModificationMotif target);
            return new Modification(_originalId: id, _modificationType: "Test", _target: target, _locationRestriction: location,
                _monoisotopicMass: mass);
        }

        internal static Protein RandomProtein(Random random, int index, int minLength, int maxLength, bool isDecoy = false)
        {
            int length = minLength + random.Next(maxLength - minLength + 1);
            char[] sequence = new char[length];
            for (int r = 0; r < length; r++)
            {
                sequence[r] = Residues[random.Next(Residues.Length)];
            }
            if (random.Next(2) == 0)
            {
                sequence[0] = 'M';
            }
            string s = new(sequence);

            // about one residue in six carries one or two modifications whose motif fits it
            var mods = new Dictionary<int, List<Modification>>();
            for (int position = 1; position <= s.Length; position++)
            {
                if (random.Next(6) != 0)
                {
                    continue;
                }
                List<Modification> fitting = Localized.Where(m => m.Target.ToString() == "X" || m.Target.ToString() == s[position - 1].ToString()).ToList();
                if (fitting.Count == 0)
                {
                    continue;
                }
                var here = new List<Modification> { fitting[random.Next(fitting.Count)] };
                Modification second = fitting[random.Next(fitting.Count)];
                if (random.Next(4) == 0 && !here.Contains(second))
                {
                    here.Add(second);
                }
                mods.Add(position, here);
            }
            return new Protein(s, (isDecoy ? "DECOY_P" : "P") + index, name: "protein " + index, isDecoy: isDecoy, oneBasedModifications: mods);
        }

        internal static DigestionParams RandomDigestionParams(Random random)
        {
            string protease = random.Next(5) switch { 0 => "top-down", 1 => "Arg-C", _ => "trypsin" };
            CleavageSpecificity mode = protease != "top-down" && random.Next(4) == 0 ? CleavageSpecificity.Semi : CleavageSpecificity.Full;
            var initiator = (InitiatorMethionineBehavior)random.Next(1, 4); // Retain, Cleave or Variable
            int isoforms = random.Next(3) == 0 ? 1 + random.Next(8) : 1024;
            return new DigestionParams(protease, random.Next(3), 1 + random.Next(8), int.MaxValue, isoforms, initiator,
                random.Next(4), mode, FragmentationTerminus.Both);
        }

        internal static (List<Modification> Fixed, List<Modification> Variable) RandomSearchMods(Random random)
        {
            var fixedMods = random.Next(3) == 0 ? new List<Modification>() : new List<Modification> { CarbamidomethylC };
            var variable = new List<Modification>();
            if (random.Next(4) != 0)
            {
                variable.Add(OxidationM);
            }
            if (random.Next(3) == 0)
            {
                variable.Add(PhosphoS);
            }
            if (random.Next(5) == 0)
            {
                variable.Add(AcetylNTerm);
            }
            return (fixedMods, variable);
        }
    }
}
