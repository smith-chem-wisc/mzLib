#nullable enable
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>
/// Builds the entrapment partner of a target protein: an isomeric sequence, its modifications
/// carried onto the residues they were on, and an accession that says what it is and where it came
/// from.
/// </summary>
public static class EntrapmentProteinGenerator
{
    /// <summary>
    /// The entrapment partner of <paramref name="target"/> for one fold.
    /// </summary>
    /// <param name="target">The target protein.</param>
    /// <param name="digestionParams">Supplies the cleavage sites, and the minimum peptide length
    /// below which a piece is not identifiable on its own.</param>
    /// <param name="forbiddenSequences">Sequences no partner peptide may equal, normally every
    /// target peptide in the database.</param>
    /// <param name="fold">Zero-based fold, in <c>[0, foldCount)</c>.</param>
    /// <param name="foldCount">Partners per target -- the <c>r</c> of an r-fold database.</param>
    /// <param name="seed">Changes every choice reproducibly.</param>
    /// <remarks>
    /// Built in two steps, both through existing constructors and neither through reflection:
    /// <see cref="Protein.CloneWithNewSequenceAndMods"/> for the sequence and modifications, then
    /// Protein's copy constructor for the accession and the entrapment flag.
    /// </remarks>
    public static Protein Create(Protein target, IDigestionParams digestionParams,
        IReadOnlySet<string> forbiddenSequences, int fold = 0, int foldCount = 1, int seed = 1,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier) =>
        Create(target, digestionParams, forbiddenSequences, out _, fold, foldCount, seed, entrapmentIdentifier);

    /// <summary>
    /// The entrapment partner of <paramref name="target"/>, along with the record of how it was
    /// built.
    /// </summary>
    /// <param name="assembly">What happened to each piece of the target: rearranged, kept because
    /// it was too short to identify, or excised for want of a partner. Feed this to
    /// <see cref="EntrapmentReportBuilder"/> rather than assembling a second time to find out.</param>
    public static Protein Create(Protein target, IDigestionParams digestionParams,
        IReadOnlySet<string> forbiddenSequences, out EntrapmentAssembly assembly,
        int fold = 0, int foldCount = 1, int seed = 1,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (target is null)
        {
            throw new MzLibException("Cannot build an entrapment protein from a null target.");
        }

        assembly = EntrapmentAssembler.Assemble(target.BaseSequence, digestionParams,
            forbiddenSequences, fold, foldCount, seed);

        Dictionary<int, List<Modification>> movedMods =
            MoveModifications(target.OneBasedPossibleLocalizedModifications, assembly.TargetToEntrapmentPosition);

        var withNewSequence = (Protein)target.CloneWithNewSequenceAndMods(assembly.EntrapmentSequence, movedMods);

        return new Protein(withNewSequence,
            accession: EntrapmentAccession.Format(target.Accession, fold, entrapmentIdentifier),
            isEntrapment: true,
            oneBasedModifications: movedMods);
    }

    /// <summary>
    /// Carries modifications onto the positions their residues moved to.
    /// </summary>
    /// <remarks>
    /// Because the rearrangement only moves positions, a modification always lands on the same
    /// residue it was on, so its motif still fits. A modification whose residue was excised has
    /// nowhere to go and is dropped rather than placed somewhere arbitrary -- silently relocating it
    /// would invent a site the target never had.
    /// </remarks>
    private static Dictionary<int, List<Modification>> MoveModifications(
        IDictionary<int, List<Modification>> oneBasedModifications, int[] targetToEntrapmentPosition)
    {
        var moved = new Dictionary<int, List<Modification>>();
        if (oneBasedModifications is null)
        {
            return moved;
        }

        foreach ((int oneBasedPosition, List<Modification> mods) in oneBasedModifications)
        {
            int zeroBased = oneBasedPosition - 1;
            if (zeroBased < 0 || zeroBased >= targetToEntrapmentPosition.Length)
            {
                continue;
            }

            int destination = targetToEntrapmentPosition[zeroBased];
            if (destination < 0)
            {
                continue;   // the residue was excised, and its modifications go with it
            }

            int newOneBased = destination + 1;
            if (moved.TryGetValue(newOneBased, out List<Modification>? existing))
            {
                existing.AddRange(mods);
            }
            else
            {
                moved[newOneBased] = new List<Modification>(mods);
            }
        }

        return moved;
    }
}
