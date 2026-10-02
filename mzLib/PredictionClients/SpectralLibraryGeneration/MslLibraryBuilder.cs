#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using MassSpectrometry;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace PredictionClients.SpectralLibraryGeneration
{
    /// <summary>How <see cref="MslLibraryBuilder"/> digests, labels and predicts.</summary>
    /// <param name="PrecursorCharges">Every peptide is predicted at each of these charges.</param>
    /// <param name="Nce">The normalized collision energy sent to the fragment model and stamped on every entry.</param>
    /// <param name="DissociationType">Sent to the fragment model as its fragmentation type and stamped on every entry.</param>
    /// <param name="DecoyType">How decoy proteins are made from the targets (and entrapment).</param>
    /// <param name="PredictionChunkSize">
    /// At most this many precursors per call to the fragment model, so a proteome-scale build is many bounded calls, with
    /// cancellation checked between them.
    /// </param>
    public sealed record MslLibraryBuildParameters(
        DigestionParams DigestionParams,
        List<Modification> FixedModifications,
        List<Modification> VariableModifications,
        IReadOnlyList<int> PrecursorCharges,
        int Nce,
        DissociationType DissociationType,
        DecoyType DecoyType = DecoyType.Reverse,
        int PredictionChunkSize = 50_000);

    /// <summary>What a build produced and what it left out, by precursor (peptide and charge).</summary>
    /// <param name="NotPredicted">Precursors a model rejected, for example too long or carrying a modification it cannot take.</param>
    /// <param name="NoFragments">Precursors predicted, but with no fragment above the intensity filter.</param>
    /// <param name="DecoysDroppedAsTargetSequences">
    /// Decoy precursors whose sequence is also a target's, with I and L counted as the same residue (same mass, same
    /// spectrum). The target keeps the sequence; such a decoy would be a target in disguise.
    /// </param>
    public sealed record MslLibraryBuildResult(int TargetPrecursors, int DecoyPrecursors, int NotPredicted, int NoFragments,
        int DecoysDroppedAsTargetSequences);

    /// <summary>
    /// Builds a predicted spectral library from proteins:
    /// <list type="number">
    /// <item>digest the targets and any entrapment proteins;</item>
    /// <item>make decoy proteins from both and digest them;</item>
    /// <item>predict fragment intensities and iRT with the given Koina models;</item>
    /// <item>return one <see cref="MslLibraryEntry"/> per peptide and charge.</item>
    /// </list>
    /// Each entry carries its decoy flag, its proteins' accessions and genes (sorted and '|'-joined when shared), the NCE and
    /// dissociation type, and the predicted iRT in its retention-time field.
    /// </summary>
    /// <remarks>
    /// Koina failures (<see cref="Koina.Client.KoinaServiceException"/>) propagate. Precursors a model rejects are counted in
    /// <see cref="MslLibraryBuildResult"/>, never thrown. Save the entries with <c>MslLibrary.Save</c>.
    /// </remarks>
    public sealed class MslLibraryBuilder
    {
        private readonly FragmentIntensityModel _intensityModel;
        private readonly RetentionTimeModel _irtModel;

        /// <exception cref="ArgumentNullException">A model is null.</exception>
        public MslLibraryBuilder(FragmentIntensityModel intensityModel, RetentionTimeModel irtModel)
        {
            _intensityModel = intensityModel ?? throw new ArgumentNullException(nameof(intensityModel));
            _irtModel = irtModel ?? throw new ArgumentNullException(nameof(irtModel));
        }

        /// <summary>
        /// Called with the target proteins before decoys are made. The proteins it returns are digested, decoyed and
        /// predicted like targets, keeping their own accessions. Null adds none.
        /// </summary>
        public Func<IReadOnlyList<Protein>, IReadOnlyList<Protein>>? EntrapmentProteins { get; init; }

        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">No precursor charge is given.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The prediction chunk size is less than 1.</exception>
        /// <exception cref="OperationCanceledException">The build was cancelled.</exception>
        public List<MslLibraryEntry> Build(IReadOnlyList<Protein> targets, MslLibraryBuildParameters parameters,
            out MslLibraryBuildResult result, CancellationToken cancellationToken = default)
        {
            ArgumentNullException.ThrowIfNull(targets);
            ArgumentNullException.ThrowIfNull(parameters);
            if (parameters.PrecursorCharges is null || parameters.PrecursorCharges.Count == 0)
                throw new ArgumentException("At least one precursor charge is needed.", nameof(parameters));
            if (parameters.PredictionChunkSize < 1)
                throw new ArgumentOutOfRangeException(nameof(parameters), parameters.PredictionChunkSize, "The prediction chunk size must be at least 1.");
            cancellationToken.ThrowIfCancellationRequested();

            var proteins = targets.ToList();
            if (EntrapmentProteins is not null)
                proteins.AddRange(EntrapmentProteins(targets));
            var decoyProteins = DecoyProteinGenerator.GenerateDecoys(proteins, parameters.DecoyType);

            var targetPeptides = Digest(proteins, parameters);
            var decoyPeptides = Digest(decoyProteins, parameters);
            // I = L: a decoy with a target's sequence up to I/L has the target's spectrum
            var targetSequences = targetPeptides.Keys.Select(LeucineEquivalent).ToHashSet(StringComparer.Ordinal);
            var clashes = decoyPeptides.Keys.Where(sequence => targetSequences.Contains(LeucineEquivalent(sequence))).ToList();
            foreach (var sequence in clashes)
                decoyPeptides.Remove(sequence);

            var precursors = new List<Precursor>();
            foreach (var (peptides, isDecoy) in new[] { (targetPeptides, false), (decoyPeptides, true) })
                foreach (var (sequence, proteinsOfPeptide) in peptides.OrderBy(p => p.Key, StringComparer.Ordinal))
                    foreach (int charge in parameters.PrecursorCharges)
                        precursors.Add(new Precursor(sequence, charge, isDecoy,
                            Join(proteinsOfPeptide.Select(p => p.Accession)),
                            Join(proteinsOfPeptide.SelectMany(p => p.GeneNames?.Select(g => g.Item2) ?? []))));

            // iRT does not depend on charge: one prediction per peptide
            var irts = PredictIrt(precursors.Select(p => p.FullSequence).Distinct(StringComparer.Ordinal).ToList(),
                parameters.PredictionChunkSize, cancellationToken);

            var entries = new List<MslLibraryEntry>();
            int notPredicted = 0, noFragments = 0;
            var predictable = new List<Precursor>();
            foreach (var precursor in precursors)
            {
                if (irts.TryGetValue(precursor.FullSequence, out double irt) && double.IsFinite(irt))
                    predictable.Add(precursor);
                else
                    notPredicted++;
            }

            string fragmentationType = parameters.DissociationType.ToString();
            foreach (var chunk in predictable.Chunk(parameters.PredictionChunkSize))
            {
                cancellationToken.ThrowIfCancellationRequested();
                _intensityModel.Predict(chunk
                    .Select(p => new FragmentIntensityPredictionInput(p.FullSequence, p.Charge, parameters.Nce, null, fragmentationType))
                    .ToList());
                notPredicted += _intensityModel.ValidInputsMask.Count(valid => !valid);

                var spectra = _intensityModel.GenerateLibrarySpectraFromPredictions(chunk.Select(p => (double?)irts[p.FullSequence]).ToArray(), out _);
                // Join by (sequence, charge): library generation de-duplicates by name, so its order is not the input's
                var byKey = chunk.ToDictionary(p => (p.FullSequence, p.Charge));
                var predicted = new HashSet<(string, int)>();
                foreach (var spectrum in spectra)
                {
                    if (!byKey.TryGetValue((spectrum.Sequence, spectrum.ChargeState), out var precursor) || !predicted.Add((spectrum.Sequence, spectrum.ChargeState)))
                        continue;
                    var entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
                    if (entry is null || entry.MatchedFragmentIons.Count == 0)
                    {
                        noFragments++;
                        continue;
                    }
                    entry.IsDecoy = precursor.IsDecoy;
                    entry.ProteinAccession = precursor.Accessions;
                    entry.GeneName = precursor.Genes;
                    entry.Nce = parameters.Nce;
                    entry.DissociationType = parameters.DissociationType;
                    entry.RetentionTime = irts[precursor.FullSequence];
                    entry.Source = MslFormat.SourceType.Predicted;
                    entries.Add(entry);
                }
            }

            result = new MslLibraryBuildResult(
                entries.Count(e => !e.IsDecoy),
                entries.Count(e => e.IsDecoy),
                notPredicted,
                noFragments,
                clashes.Count * parameters.PrecursorCharges.Count);
            return entries;
        }

        /// <summary>The full sequence with every residue I read as L; modification names in brackets are left alone.</summary>
        private static string LeucineEquivalent(string fullSequence)
        {
            var chars = fullSequence.ToCharArray();
            int depth = 0;
            for (int i = 0; i < chars.Length; i++)
            {
                if (chars[i] == '[') depth++;
                else if (chars[i] == ']') depth--;
                else if (depth == 0 && chars[i] == 'I') chars[i] = 'L';
            }
            return new string(chars);
        }

        private sealed record Precursor(string FullSequence, int Charge, bool IsDecoy, string Accessions, string Genes);

        /// <summary>Each digested peptide (by full sequence, finite mass only) with every protein it comes from.</summary>
        private static Dictionary<string, List<Protein>> Digest(IEnumerable<Protein> proteins, MslLibraryBuildParameters parameters)
        {
            var peptides = new Dictionary<string, List<Protein>>(StringComparer.Ordinal);
            foreach (var protein in proteins)
                foreach (var peptide in protein.Digest(parameters.DigestionParams, parameters.FixedModifications, parameters.VariableModifications))
                {
                    if (!double.IsFinite(peptide.MonoisotopicMass))
                        continue;
                    if (!peptides.TryGetValue(peptide.FullSequence, out var list))
                        peptides[peptide.FullSequence] = list = new List<Protein>();
                    if (!list.Contains(protein))
                        list.Add(protein);
                }
            return peptides;
        }

        private Dictionary<string, double> PredictIrt(List<string> sequences, int chunkSize, CancellationToken cancellationToken)
        {
            var irts = new Dictionary<string, double>(StringComparer.Ordinal);
            foreach (var chunk in sequences.Chunk(chunkSize))
            {
                cancellationToken.ThrowIfCancellationRequested();
                var predictions = _irtModel.Predict(chunk.Select(s => new RetentionTimePredictionInput(s)).ToList());
                for (int i = 0; i < chunk.Length; i++)
                    if (predictions[i].PredictedRetentionTime is double irt)
                        irts[chunk[i]] = irt;
            }
            return irts;
        }

        private static string Join(IEnumerable<string> values) =>
            string.Join("|", values.Where(v => !string.IsNullOrEmpty(v)).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal));
    }
}
