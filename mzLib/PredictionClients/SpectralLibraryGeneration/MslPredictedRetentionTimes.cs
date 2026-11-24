#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using Chromatography.RetentionTimeCalibration;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using MzLibUtil;
using Omics.SequenceConversion;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using Readers.ProForma;

namespace PredictionClients.SpectralLibraryGeneration
{
    /// <summary>Which entries <see cref="MslPredictedRetentionTimes.Assign"/> gives a predicted iRT.</summary>
    public enum MslPredictedIrtMode
    {
        /// <summary>
        /// Every entry. For a library whose stored times cannot be used, such as one converted from
        /// <c>.msp</c>: an <c>.msp</c> time is minutes in one unnamed run, which no calibration can map
        /// onto the iRT scale. Every status is reset to none.
        /// </summary>
        ReplaceAll,

        /// <summary>
        /// Only entries with no time (NaN) and no status (never measured). Never touches a measured time,
        /// so it is safe to repeat, for example after a predictor outage left some times missing.
        /// </summary>
        FillMissing
    }

    /// <summary>What <see cref="MslPredictedRetentionTimes.Assign"/> did.</summary>
    /// <param name="Considered">Entries the mode selected.</param>
    /// <param name="FromIrtPredictor">Of those, times from the iRT predictor (peptides of up to 30 residues).</param>
    /// <param name="FromLongPeptidePredictor">Times from the long-peptide predictor (unmodified, 31-50 residues).</param>
    /// <param name="NoPrediction">Entries left with no time (NaN): too long, a modification the predictor refuses, or no predictor.</param>
    /// <param name="SequencesPredicted">Distinct full sequences sent to the predictors (iRT does not depend on charge).</param>
    public sealed record MslPredictedIrtResult(int Considered, int FromIrtPredictor, int FromLongPeptidePredictor, int NoPrediction,
        int SequencesPredicted);

    /// <summary>
    /// Gives library entries a predicted iRT: the start of a calibrated library. A predicted time has status
    /// none, so the first measured time replaces it (<c>MslLibrary.UpdateRetentionTimesAndSave</c>).
    /// </summary>
    /// <example>
    /// Converting a MetaMorpheus <c>.msp</c> library, whose times are run minutes:
    /// <code>
    /// var prosit = MslPredictedRetentionTimes.FromKoina(new Prosit2025iRT40PTM());
    /// List&lt;MslLibraryEntry&gt; entries = MslConverter.FromMspFile(mspPath).GetAllEntries().ToList();
    /// entries.ForEach(e =&gt; e.Source = MslFormat.SourceType.Empirical);   // the spectra are measured
    /// entries = MslPredictedRetentionTimes.Assign(entries, prosit, MslPredictedIrtMode.ReplaceAll, out _);
    ///
    /// // Optionally, unmodified peptides of 31-50 residues from Chronologer, mapped onto Prosit's scale
    /// var known = entries.Where(e =&gt; double.IsFinite(e.RetentionTime))
    ///                    .GroupBy(e =&gt; e.FullSequence).ToDictionary(g =&gt; g.Key, g =&gt; g.First().RetentionTime);
    /// var chronologer = MslPredictedRetentionTimes.ChronologerOnIrtScale(known);
    /// entries = MslPredictedRetentionTimes.Assign(entries, prosit, MslPredictedIrtMode.FillMissing, out _, chronologer);
    /// MslLibrary.Save(mslPath, entries);
    /// </code>
    /// </example>
    public static class MslPredictedRetentionTimes
    {
        private const string FixedCarbamidomethyl = "Common Fixed:Carbamidomethyl on C";

        /// <summary>A bracketed modification label, allowing one level of nested brackets.</summary>
        private static readonly Regex ModificationLabel = new(@"\[((?:[^\[\]]|\[[^\]]*\])*)\]", RegexOptions.Compiled);

        /// <summary>The longest peptide the iRT predictor is asked for (Prosit's limit).</summary>
        public const int MaximumIrtPredictorLength = 30;

        /// <summary>The longest peptide the long-peptide predictor is asked for (Chronologer's limit).</summary>
        public const int MaximumLongPeptideLength = 50;

        /// <summary>
        /// Returns copies of <paramref name="entries"/> in which the entries <paramref name="mode"/> selects
        /// have a predicted iRT; the others are returned as they are, and no input entry is modified.
        ///
        /// <para>
        /// A peptide of up to 30 residues (any modifications) gets <paramref name="predictIrt"/>'s value; an
        /// unmodified one of 31-50 residues gets <paramref name="predictLongPeptideIrt"/>'s, when given; any
        /// other gets NaN. Fixed carbamidomethyl C counts as unmodified. Each distinct full sequence is
        /// predicted once, whatever its charges.
        /// </para>
        ///
        /// <para>
        /// A changed entry has <see cref="MslLibraryEntry.RtStatus"/> none, no strike, and
        /// <see cref="MslLibraryEntry.RtIsCalibrated"/> false (the time is iRT). Its source is unchanged: set
        /// <see cref="MslFormat.SourceType.Empirical"/> yourself when the spectra are measured.
        /// </para>
        /// </summary>
        /// <param name="entries">The library's entries.</param>
        /// <param name="predictIrt">
        ///   iRT for each full sequence given, in order; NaN where there is none. See <see cref="FromKoina"/>.
        /// </param>
        /// <param name="mode">Which entries to give a time.</param>
        /// <param name="result">What was done.</param>
        /// <param name="predictLongPeptideIrt">
        ///   Optional iRT for unmodified peptides of 31-50 residues. See <see cref="ChronologerOnIrtScale"/>.
        /// </param>
        /// <exception cref="ArgumentNullException"><paramref name="entries"/> or <paramref name="predictIrt"/> is null.</exception>
        /// <exception cref="ArgumentException">An entry, or its sequence, is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="mode"/> is not defined.</exception>
        /// <exception cref="InvalidOperationException">A predictor returned null or the wrong number of values.</exception>
        public static List<MslLibraryEntry> Assign(IReadOnlyList<MslLibraryEntry> entries,
            Func<IReadOnlyList<string>, IReadOnlyList<double>> predictIrt, MslPredictedIrtMode mode,
            out MslPredictedIrtResult result,
            Func<IReadOnlyList<string>, IReadOnlyList<double>>? predictLongPeptideIrt = null)
        {
            ArgumentNullException.ThrowIfNull(entries);
            ArgumentNullException.ThrowIfNull(predictIrt);
            if (!Enum.IsDefined(mode))
                throw new ArgumentOutOfRangeException(nameof(mode), mode, "Not a defined MslPredictedIrtMode.");
            for (int i = 0; i < entries.Count; i++)
                if (entries[i] is null || entries[i].FullSequence is null || entries[i].BaseSequence is null)
                    throw new ArgumentException($"entries[{i}] is null or has a null sequence.", nameof(entries));

            // The entries to change, routed to a predictor by length and modification
            var selected = new List<int>();
            var forIrt = new HashSet<string>(StringComparer.Ordinal);
            var forLong = new HashSet<string>(StringComparer.Ordinal);
            for (int i = 0; i < entries.Count; i++)
            {
                MslLibraryEntry e = entries[i];
                bool select = mode == MslPredictedIrtMode.ReplaceAll ||
                    (e.RtStatus == MslFormat.RetentionTimeStatus.None && double.IsNaN(e.RetentionTime));
                if (!select)
                    continue;
                selected.Add(i);
                int length = e.BaseSequence.Length;
                if (length <= MaximumIrtPredictorLength)
                    forIrt.Add(e.FullSequence);
                else if (predictLongPeptideIrt is not null && length <= MaximumLongPeptideLength && IsUnmodified(e.FullSequence))
                    forLong.Add(e.FullSequence);
            }

            Dictionary<string, double> irt = Predict(predictIrt, forIrt);
            Dictionary<string, double> longIrt = predictLongPeptideIrt is null
                ? new Dictionary<string, double>()
                : Predict(predictLongPeptideIrt, forLong);

            var output = entries.ToList();
            int fromIrt = 0, fromLong = 0, none = 0;
            foreach (int i in selected)
            {
                MslLibraryEntry e = entries[i];
                double time = double.NaN;
                if (irt.TryGetValue(e.FullSequence, out double p) && double.IsFinite(p)) { time = p; fromIrt++; }
                else if (longIrt.TryGetValue(e.FullSequence, out double q) && double.IsFinite(q)) { time = q; fromLong++; }
                else none++;

                MslLibraryEntry copy = e.WithFragments(e.MatchedFragmentIons);
                copy.RetentionTime = time;
                copy.RtStatus = MslFormat.RetentionTimeStatus.None;
                copy.RtStruckByPrediction = false;
                copy.RtIsCalibrated = false;
                output[i] = copy;
            }

            result = new MslPredictedIrtResult(selected.Count, fromIrt, fromLong, none, forIrt.Count + forLong.Count);
            return output;
        }

        /// <summary>
        /// Wraps a Koina iRT model (for example <c>Prosit2025iRT40PTM</c>) as a predictor for
        /// <see cref="Assign"/> and for <c>MslLibrary.UpdateRetentionTimesAndSave</c>.
        ///
        /// <para>
        /// Each full sequence (mzLib notation) is converted to ProForma, so a catalog modification reaches the
        /// model by its Unimod id, whatever name it carries. A sequence that cannot be converted, or that the
        /// model refuses (a residue or modification it cannot take, or more than 30 residues), gets NaN.
        /// </para>
        /// </summary>
        /// <param name="model">An iRT model (<see cref="RetentionTimeModel.IsIndexedRetentionTimeModel"/>).</param>
        /// <param name="chunkSize">At most this many sequences per call to the model; cancellation is checked between calls.</param>
        /// <param name="cancellationToken">Stops between chunks.</param>
        /// <returns>A predictor; an outage of the service propagates from it as the model throws it.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="model"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="model"/> predicts minutes, not iRT.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="chunkSize"/> is less than 1.</exception>
        public static Func<IReadOnlyList<string>, IReadOnlyList<double>> FromKoina(RetentionTimeModel model,
            int chunkSize = 50_000, CancellationToken cancellationToken = default)
        {
            ArgumentNullException.ThrowIfNull(model);
            if (!model.IsIndexedRetentionTimeModel)
                throw new ArgumentException("The model must predict iRT, not run-specific minutes.", nameof(model));
            ArgumentOutOfRangeException.ThrowIfLessThan(chunkSize, 1);
            ProFormaSequenceConversion.RegisterWithDefault();

            return sequences =>
            {
                var values = Enumerable.Repeat(double.NaN, sequences.Count).ToArray();
                var inputs = new List<(int Index, string ProForma)>();
                for (int i = 0; i < sequences.Count; i++)
                {
                    string? proForma = SequenceConversionService.Default.Convert(sequences[i], "mzLib", "ProForma",
                        null, SequenceConversionHandlingMode.ReturnNull);
                    if (proForma is not null)
                        inputs.Add((i, proForma));
                }

                foreach (var chunk in inputs.Chunk(chunkSize))
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    var predictions = model.Predict(chunk
                        .Select(x => new RetentionTimePredictionInput(x.ProForma) { SequenceParser = ProFormaSequenceParser.Instance })
                        .ToList());
                    for (int k = 0; k < chunk.Length; k++)
                        if (predictions[k].PredictedRetentionTime is double v && double.IsFinite(v))
                            values[chunk[k].Index] = v;
                }
                return values;
            };
        }

        /// <summary>
        /// A predictor for unmodified peptides of 31-50 residues, beyond Prosit's limit: Chronologer (run
        /// locally) predicts % acetonitrile, which is mapped onto the iRT scale by a robust fit
        /// (<see cref="IrtCalibration"/>, LOWESS) over peptides whose iRT is known.
        /// </summary>
        /// <param name="knownIrt">
        ///   Full sequences with a known iRT, for example a library's Prosit times. Only unmodified ones are
        ///   used, at most <paramref name="maximumFitPeptides"/> of them, chosen by SHA-256 of the sequence so
        ///   the fit does not depend on input order.
        /// </param>
        /// <param name="maximumFitPeptides">The most peptides Chronologer is run on for the fit.</param>
        /// <returns>A predictor: the mapped iRT of each unmodified sequence of up to 50 residues, NaN otherwise.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="knownIrt"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="maximumFitPeptides"/> is less than 20.</exception>
        /// <exception cref="ArgumentException">Too few known peptides for a fit (IrtCalibration needs 20).</exception>
        public static Func<IReadOnlyList<string>, IReadOnlyList<double>> ChronologerOnIrtScale(
            IReadOnlyDictionary<string, double> knownIrt, int maximumFitPeptides = 10_000)
        {
            ArgumentNullException.ThrowIfNull(knownIrt);
            ArgumentOutOfRangeException.ThrowIfLessThan(maximumFitPeptides, 20);

            var fitSet = knownIrt
                .Where(kv => kv.Key is not null && double.IsFinite(kv.Value) && IsUnmodified(kv.Key) && BaseOf(kv.Key).Length <= MaximumLongPeptideLength)
                .OrderBy(kv => Sha256Hex(kv.Key), StringComparer.Ordinal)
                .Take(maximumFitPeptides)
                .ToList();
            double[] acn = Chronologer(fitSet.Select(kv => kv.Key).ToList());
            var anchors = new List<(RtMinutes, Irt)>();
            for (int i = 0; i < fitSet.Count; i++)
                if (double.IsFinite(acn[i]))
                    anchors.Add((new RtMinutes(acn[i]), new Irt(fitSet[i].Value))); // % ACN in the "minutes" slot

            IrtCalibrationModel map = IrtCalibration.Fit(anchors, new IrtCalibrationOptions(MinimumBinShare: 0));

            return sequences =>
            {
                double[] values = Chronologer(sequences);
                for (int i = 0; i < values.Length; i++)
                    values[i] = double.IsFinite(values[i]) ? map.ToIrt(new RtMinutes(values[i])).Value : double.NaN;
                return values;
            };
        }

        /// <summary>Chronologer % ACN for each unmodified sequence of up to 50 residues; NaN for any other.</summary>
        private static double[] Chronologer(IReadOnlyList<string> sequences)
        {
            var values = Enumerable.Repeat(double.NaN, sequences.Count).ToArray();
            var asked = new Dictionary<string, IRetentionPredictable>(StringComparer.Ordinal);
            foreach (string s in sequences)
            {
                if (s is null || asked.ContainsKey(s) || !IsUnmodified(s)) continue;
                string bare = BaseOf(s);
                if (bare.Length <= MaximumLongPeptideLength)
                    asked[s] = new SequenceOnly(s, bare);
            }
            if (asked.Count == 0)
                return values;

            var acn = new Dictionary<string, double>(StringComparer.Ordinal);
            using (var chronologer = new ChronologerRetentionTimePredictor(SequenceConversionHandlingMode.ReturnNull))
                foreach (var (value, peptide, _) in chronologer.PredictRetentionTimeEquivalents(asked.Values.ToList()))
                    if (value is double v && double.IsFinite(v))
                        acn[peptide.FullSequence] = v;

            for (int i = 0; i < sequences.Count; i++)
                if (sequences[i] is not null && acn.TryGetValue(sequences[i], out double v))
                    values[i] = v;
            return values;
        }

        /// <summary>A peptide as the predictors need it: only its sequences.</summary>
        private sealed record SequenceOnly(string FullSequence, string BaseSequence) : IRetentionPredictable
        {
            public double MonoisotopicMass => 0;
            public string FullSequenceWithMassShifts => FullSequence;
        }

        private static Dictionary<string, double> Predict(Func<IReadOnlyList<string>, IReadOnlyList<double>> predict,
            HashSet<string> sequences)
        {
            var values = new Dictionary<string, double>(StringComparer.Ordinal);
            if (sequences.Count == 0)
                return values;
            var asked = sequences.Order(StringComparer.Ordinal).ToList();
            IReadOnlyList<double> got = predict(asked)
                ?? throw new InvalidOperationException("The retention-time predictor returned null.");
            if (got.Count != asked.Count)
                throw new InvalidOperationException(
                    $"The retention-time predictor returned {got.Count} values for {asked.Count} sequences.");
            for (int i = 0; i < asked.Count; i++)
                values[asked[i]] = double.IsFinite(got[i]) ? got[i] : double.NaN;
            return values;
        }

        /// <summary>No modification but fixed carbamidomethyl C.</summary>
        private static bool IsUnmodified(string fullSequence) =>
            ModificationLabel.Matches(fullSequence).All(m => m.Groups[1].Value == FixedCarbamidomethyl);

        private static string BaseOf(string fullSequence) => ModificationLabel.Replace(fullSequence, "");

        private static string Sha256Hex(string s) =>
            Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(s))).ToLowerInvariant();
    }
}
