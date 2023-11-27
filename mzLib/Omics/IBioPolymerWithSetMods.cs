using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
    public interface IBioPolymerWithSetMods : IHasChemicalFormula
    {
        string BaseSequence { get; }
        string FullSequence { get; }
        double MostAbundantMonoisotopicMass { get; }
        string SequenceWithChemicalFormulas { get; }

        int OneBasedStartResidue { get; }
        int OneBasedEndResidue { get; }
        int MissedCleavages { get; }
        CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; }
        char PreviousResidue { get; }
        char NextResidue { get; }



        //IDigestionParams DigestionParams { get; }



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

        public string DetermineFullSequence()
        {
            var subSequence = new StringBuilder();

            // modification on peptide N-terminus
            if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
            {
                subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subSequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                subSequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            return subSequence.ToString();
        }
    }
}
