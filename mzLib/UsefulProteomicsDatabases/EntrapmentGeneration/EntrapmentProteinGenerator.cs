#nullable enable
using System;
using MzLibUtil;
using Omics.Digestion;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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
    /// Entrapment entries taken from a foreign proteome, and the peptides they share with the
    /// target database.
    /// </summary>
    /// <param name="foreignProteins">Proteins of a species absent from the sample.</param>
    /// <param name="digestionParams">Used to digest the foreign proteins for the sharing check.</param>
    /// <param name="targetPeptides">Every peptide of the target database.</param>
    /// <param name="sharedWithTarget">Foreign peptides that are also target peptides, by the
    /// foreign protein's own accession. Empty when the two proteomes are disjoint.</param>
    /// <remarks>
    /// <para><b>Why this arm is required scope rather than a nicety.</b> If the entrapment sequences
    /// and the decoy sequences are both shuffles of the target, they are the same construction, and
    /// comparing one against the other demonstrates nothing about either -- the circularity Madej
    /// and Lam describe. A foreign proteome is the one entrapment source that is not derived from
    /// the target at all, so it is what makes that comparison mean something.</para>
    /// <para><b>Nothing is rearranged.</b> A foreign protein already is a sequence the sample cannot
    /// contain, so it is relabelled, not permuted; its own annotations describe its own sequence and
    /// stay valid, unlike the permutation path where every positional annotation had to be dropped.
    /// It has no target partner by construction, and <see cref="EntrapmentAccession.TryParse"/>
    /// refuses its accession for that reason.</para>
    /// <para><b>The hazard is homology, and it is counted rather than repaired.</b> A conserved
    /// protein shares peptides with its ortholog, and a shared peptide is a <i>real target peptide</i>
    /// sitting in the entrapment database -- a search matching it counts a true peptide as an
    /// entrapment discovery. Nothing can be permuted away here: the sequence is what it is. So they
    /// are reported, and the caller decides whether to exclude the peptides or drop the proteins.
    /// That is the same rule this generator applies everywhere: flag and count, never silently
    /// repair.</para>
    /// </remarks>
    public static List<Protein> GenerateForeignEntrapment(IEnumerable<Protein> foreignProteins,
        IDigestionParams digestionParams, IReadOnlySet<string> targetPeptides,
        out IReadOnlyDictionary<string, IReadOnlyCollection<string>> sharedWithTarget,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (foreignProteins is null)
        {
            throw new MzLibException("Cannot build a foreign-species arm from a null protein list.");
        }
        if (targetPeptides is null)
        {
            throw new MzLibException("The sharing check needs the target database's peptides.");
        }

        var shared = new Dictionary<string, IReadOnlyCollection<string>>();
        var entrapment = new List<Protein>();
        var noMods = new List<Modification>();

        foreach (Protein foreign in DatabaseEntries(foreignProteins))
        {
            var collisions = new HashSet<string>();
            foreach (var peptide in foreign.Digest(digestionParams, noMods, noMods))
            {
                if (targetPeptides.Contains(peptide.BaseSequence))
                {
                    collisions.Add(peptide.BaseSequence);
                }
            }

            if (collisions.Count > 0)
            {
                shared[foreign.Accession] = collisions;
            }

            entrapment.Add(CreateForeign(foreign, entrapmentIdentifier));
        }

        sharedWithTarget = shared;
        return entrapment;
    }

    /// <summary>One foreign protein, relabelled as an entrapment entry.</summary>
    public static Protein CreateForeign(Protein foreign,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (foreign is null)
        {
            throw new MzLibException("Cannot build a foreign entrapment entry from a null protein.");
        }

        // The sequence is untouched, so unlike the permutation path the positional annotations still
        // describe it and are kept. Sequence variations are the exception: applying them would
        // expand one entry into several, and the entrapment side is one entry per entry.
        return new Protein(foreign,
            accession: EntrapmentAccession.FormatForeign(foreign.Accession, entrapmentIdentifier),
            isEntrapment: true,
            sequenceVariations: new List<SequenceVariation>(),
            appliedSequenceVariations: new List<SequenceVariation>());
    }

    /// <summary>
    /// The entrapment partners for a whole target database: one per database <b>entry</b>, per fold.
    /// </summary>
    /// <param name="targets">The target proteins, as a loader returns them.</param>
    /// <param name="digestionParams">Supplies the cleavage sites, and the minimum peptide length
    /// below which a piece is not identifiable on its own.</param>
    /// <param name="forbiddenSequences">Sequences no partner peptide may equal, normally every
    /// target peptide in the database.</param>
    /// <param name="foldCount">Partners per target -- the <c>r</c> of an r-fold database.</param>
    /// <param name="seed">Changes every choice reproducibly.</param>
    /// <remarks>
    /// <para><b>Why this exists, and why it is not a <c>foreach</c> over the list.</b> A protein
    /// database is a list of entries, but a loader does not hand back a list of entries: it applies
    /// each entry's sequence variations and hands back one protein per applied combination, all of
    /// them sharing one <see cref="Protein.ConsensusVariant"/>. Loading the reviewed human XML gives
    /// 52,337 proteins from 20,416 entries.</para>
    /// <para>Generating a partner for each of those proteins produces 52,337 entrapment proteins,
    /// and because each is constructed fresh it is its own consensus.
    /// <see cref="ProteinDbWriter.WriteXmlDatabase"/> writes one entry per distinct consensus, so
    /// the targets fold back to 20,416 entries and the entrapment proteins do not fold at all. The
    /// database that reaches a search engine then holds <b>2.56 entrapment entries per target</b>,
    /// while every FDP estimator reading it is told <c>r = 1</c>. Measured, not hypothetical: it is
    /// what the human XML arm of this project's own experiment was searched against.</para>
    /// <para>So the unit of entrapment is the entry. A sequence variant is an annotation on an
    /// entry, not an entry of its own, and an entrapment partner for every point variant of a
    /// protein is not something anyone asked for -- it is only what a naive loop produces. The
    /// FASTA path is unaffected either way, because a FASTA has no variants and every protein is
    /// already its own consensus.</para>
    /// </remarks>
    public static List<Protein> GenerateEntrapment(IEnumerable<Protein> targets,
        IDigestionParams digestionParams, IReadOnlySet<string> forbiddenSequences,
        int foldCount = 1, int seed = 1,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
        => GenerateEntrapment(targets, digestionParams, forbiddenSequences, null, foldCount, seed,
            entrapmentIdentifier);

    /// <summary>
    /// As above, reporting each partner's <see cref="EntrapmentAssembly"/> as it is produced.
    /// </summary>
    /// <param name="observe">Called with the entry, the fold, and how that partner was built.</param>
    /// <remarks>
    /// The parameterless overload discards the assembly, so a caller building a QC report had to
    /// write its own <see cref="DatabaseEntries"/> plus per-protein <see cref="Create(Protein, IDigestionParams, IReadOnlySet{string}, out EntrapmentAssembly, int, int, int, string)"/>
    /// loop -- two loops obliged to agree about folds and about what an entry is. They stopped
    /// agreeing: one such caller compared its partner count against the <i>loaded</i> protein count
    /// rather than the entry count and failed silently for a day. One loop cannot disagree with
    /// itself.
    /// </remarks>
    public static List<Protein> GenerateEntrapment(IEnumerable<Protein> targets,
        IDigestionParams digestionParams, IReadOnlySet<string> forbiddenSequences,
        Action<Protein, int, EntrapmentAssembly>? observe,
        int foldCount = 1, int seed = 1,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (targets is null)
        {
            throw new MzLibException("Cannot build entrapment proteins from a null target list.");
        }

        // Refused here as well as inside the assembler. The assembler's copy fires on the first
        // protein so nothing is wasted either way, but a throw from that depth reaches a GUI user
        // as a crash part-way through a run rather than as "this agent cannot be used", and this
        // is the call they actually made.
        EntrapmentAssembler.RefuseAgentsWhoseSitesCannotBeHeld(
            digestionParams?.DigestionAgent?.Name ?? string.Empty,
            digestionParams?.DigestionAgent?.DigestionMotifs ?? new List<DigestionMotif>());

        var entrapment = new List<Protein>();
        foreach (Protein entry in DatabaseEntries(targets))
        {
            for (int fold = 0; fold < foldCount; fold++)
            {
                Protein partner = Create(entry, digestionParams, forbiddenSequences, out EntrapmentAssembly assembly,
                    fold, foldCount, seed, entrapmentIdentifier);

                // Reported before the emptiness check, so a report still sees a protein that
                // produced nothing. An entry that contributes no partner is exactly what a reader of
                // the achieved ratio needs to know about.
                observe?.Invoke(entry, fold, assembly);

                // Every piece excised leaves nothing behind. A protein with an empty sequence is not
                // a database entry: a loader warns and discards it, and its decoy with it.
                // Q156A1 is the real case -- a methionine followed by 79 glutamines, whose only base
                // piece has exactly one arrangement once the termini are anchored.
                //
                // The loaded database is unchanged by skipping it: the loader discarded it anyway,
                // and the target, entrapment and decoy counts come out the same either way. What
                // skipping it removes is the loader's "N empty entries ignored" warning, which was
                // the only signal that the arm was short. EntrapmentReport.EntriesYieldingNoPartner
                // carries that instead -- short by construction, said where it is decided, rather
                // than short and then noticed by whatever happens to load the file.
                if (assembly.EntrapmentSequence.Length == 0)
                {
                    continue;
                }

                entrapment.Add(partner);
            }
        }

        return entrapment;
    }

    /// <summary>
    /// The entries behind a loaded protein list: each distinct consensus, once, in the order first
    /// seen. This is the list an entrapment database is one-to-one with, and the list a database
    /// writer will persist.
    /// </summary>
    /// <remarks>
    /// <para>Public because <see cref="GenerateEntrapment"/> cannot serve every caller: building the
    /// QC report needs the <see cref="EntrapmentAssembly"/> of each partner, which only the
    /// per-protein <see cref="Create(Protein, IDigestionParams, IReadOnlySet{string}, out EntrapmentAssembly, int, int, int, string)"/>
    /// overload returns. Such a caller must iterate the same entries, and should not have to
    /// re-derive what an entry is.</para>
    /// <para>Keyed by reference, not by value. <see cref="Protein"/> overrides <c>Equals</c>, and the
    /// question here is which entry a protein came from -- object identity -- rather than whether
    /// two proteins happen to look alike. <see cref="DecoyProteinGenerator"/> keys its own consensus
    /// map the same way and for the same reason.</para>
    /// </remarks>
    public static IEnumerable<Protein> DatabaseEntries(IEnumerable<Protein> targets)
    {
        var seen = new HashSet<object>(ReferenceEqualityComparer.Instance);
        foreach (Protein target in targets)
        {
            if (target is null)
            {
                continue;
            }

            // A protein with no variants is its own consensus, so this is the protein itself.
            Protein entry = target.ConsensusVariant as Protein ?? target;
            if (seen.Add(entry))
            {
                yield return entry;
            }
        }
    }

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
            uniProtSequenceAttributes: DescribeSequence(target, assembly.EntrapmentSequence),
            oneBasedModifications: movedMods,
            // Every other positional annotation describes the TARGET's sequence and means nothing
            // once the residues have moved. Carrying them over is not merely untidy: excision
            // shortens the protein, so a coordinate can point past its end, and a consumer that
            // indexes with it -- MetaMorpheus applies sequence variations while loading a database
            // -- throws. Measured on the human proteome: 131 entrapment entries carried a feature
            // position beyond their own length, and every search against that database failed.
            // Modifications are the one annotation that survives, because they alone are remapped.
            sequenceVariations: new List<SequenceVariation>(),
            appliedSequenceVariations: new List<SequenceVariation>(),
            proteolysisProducts: new List<TruncationProduct>(),
            disulfideBonds: new List<DisulfideBond>(),
            spliceSites: new List<SpliceSite>());
    }

    /// <summary>
    /// Sequence attributes describing the entrapment sequence rather than the one it came from.
    /// </summary>
    /// <remarks>
    /// The copy constructor carries the target's <c>UniProtSequenceAttributes</c> across, and
    /// excision changes the length -- so a 15-residue partner was written as
    /// <c>&lt;sequence length="30" mass="3114"&gt;</c>. <see cref="ProteinDbLoader.LoadProteinXML"/>
    /// recomputes length on read, which is why a round trip through mzLib never noticed, but the file
    /// on disk was self-inconsistent for anything else that reads it.
    /// <para>Mass is recomputed too, though for a different reason: a rearrangement is isomeric, so
    /// an unexcised partner has exactly its target's mass and only an excised one does not. Deriving
    /// both from the sequence is right in either case and does not depend on knowing which.</para>
    /// </remarks>
    private static UniProtSequenceAttributes? DescribeSequence(Protein target, string entrapmentSequence)
    {
        UniProtSequenceAttributes? source = target.UniProtSequenceAttributes;
        if (source is null)
        {
            return null;   // the target carried none, so inventing some would assert more than is known
        }

        // Checksum and sequence version describe the TARGET's sequence, and the partner's is a
        // different one. Carrying them across wrote every entrapment entry with its target's CRC64
        // -- a checksum that both ProteinXmlEntry parses and ProteinDbWriter writes, so it is read
        // back as an assertion about a sequence it does not describe. "unknown" and an unset
        // version say the true thing, and are what Protein's own constructor synthesizes when it is
        // handed nothing.
        var described = new UniProtSequenceAttributes(source.Length, source.Mass,
            checkSum: "unknown", source.EntryModified, sequenceVersion: -1,
            source.IsPrecursor, source.Fragment);
        described.UpdateLengthAttribute(entrapmentSequence);

        // MonoisotopicMass is NaN for a sequence holding X or B, and (int)Math.Round(NaN) is 0, so
        // an entry whose target carried a real mass was written as mass="0". A rearrangement is
        // isomeric, so the target's own mass is the right answer whenever the partner was not
        // excised; where it was, an inherited mass is still closer to true than zero.
        double mass = new PeptideWithSetModifications(entrapmentSequence,
            new Dictionary<string, Modification>()).MonoisotopicMass;
        if (double.IsFinite(mass))
        {
            described.UpdateMassAttribute(entrapmentSequence);
        }

        return described;
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
