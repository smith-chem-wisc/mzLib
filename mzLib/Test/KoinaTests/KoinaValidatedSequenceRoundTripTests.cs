using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// ValidatedFullSequence and the Koina payload must name the same modifications. Every protein modification in
    /// mzLib's catalogs is put on a peptide by digestion, wherever its motif and location restriction let it go, and the
    /// peptide is written by mzLib's writers (full sequence, ProForma, mass shifts). On one model per distinct
    /// modification policy: an accepted ValidatedFullSequence is the full sequence or ProForma string as written (its
    /// modifications keep their text), re-cleaning it with the same parser gives the same payload, and the peptide
    /// built from it, and from the input, is the digested peptide: its mass, its modifications at the same places, and
    /// for the full sequence its label.
    /// </summary>
    [TestFixture]
    public class KoinaValidatedSequenceRoundTripTests
    {
        private const string Residues = "ACDEFGHIKLMNPQRSTVWY";
        private const double MassTolerance = 0.01;

        private static readonly Lazy<List<(PeptideWithSetModifications Peptide, double ModificationMass)>> Digests = new(() =>
            Mods.AllProteinModsList
                .Where(m => m.MonoisotopicMass.HasValue && m.Target != null && m.Target.ToString().Length == 1
                    && (Residues.Contains(m.Target.ToString()) || m.Target.ToString() == "X"))
                .SelectMany(m =>
                {
                    var target = m.Target.ToString() == "X" ? "A" : m.Target.ToString();
                    var protein = m.LocationRestriction switch
                    {
                        "N-terminal." => $"{target}GLSDEFK",
                        "C-terminal." => $"GLSDEF{target}",
                        "Peptide C-terminal." => $"GLSDEF{target}AGLSDEFK",
                        _ => $"GLS{target}DEFK"
                    };
                    return new Protein(protein, "P")
                        .Digest(new DigestionParams("trypsin", maxMissedCleavages: 1, minPeptideLength: 1, maxModsForPeptides: 1),
                            new List<Modification>(), new List<Modification> { m })
                        .Where(p => p.AllModsOneIsNterminus.Count == 1)
                        .Select(p => (p, m.MonoisotopicMass!.Value));
                })
                .DistinctBy(d => d.p.FullSequence)
                .ToList());

        [Test]
        public void ValidatedFullSequence_RoundTripsToTheDigestedPeptide_FullSequence() =>
            AssertRoundTrips(p => p.FullSequence, null, minimumAccepted: 1000, asWritten: true);

        [Test]
        public void ValidatedFullSequence_RoundTripsToTheDigestedPeptide_ProForma() =>
            AssertRoundTrips(p => p.ToProFormaString(), ProFormaSequenceParser.Instance, minimumAccepted: 1000, validatedAsWritten: true);

        [Test]
        public void ValidatedFullSequence_RoundTripsToTheDigestedPeptide_MassShifts() =>
            AssertRoundTrips(p => p.FullSequenceWithMassShifts, MassShiftSequenceParser.Instance, minimumAccepted: 1000);

        // Validation is the same in every model; what differs is the modification policy, so one model per policy.
        // Fragment models carry every policy except detectability's, which allows no modification, and they alone build
        // peptides. A model with an allow-list is only given the modifications whose mass one of its ids could match.
        private static IEnumerable<(object Model, Func<double, bool> CouldAccept)> ModelsByPolicy()
        {
            var unimodMasses = Mods.UnimodModifications
                .Where(m => m.MonoisotopicMass.HasValue && m.DatabaseReference != null)
                .SelectMany(m => m.DatabaseReference
                    .Where(reference => reference.Key.Equals("Unimod", StringComparison.OrdinalIgnoreCase))
                    .SelectMany(reference => reference.Value)
                    .Select(id => (Ok: int.TryParse(id, out var unimodId), Id: unimodId, Mass: m.MonoisotopicMass!.Value)))
                .Where(entry => entry.Ok)
                .GroupBy(entry => entry.Id)
                .ToDictionary(g => g.Key, g => g.First().Mass);

            var models = typeof(FragmentIntensityModel).Assembly.GetTypes()
                .Where(t => !t.IsAbstract && typeof(FragmentIntensityModel).IsAssignableFrom(t))
                .OrderBy(t => t.FullName)
                .Select(t =>
                {
                    var ctor = t.GetConstructors().First(c => c.GetParameters().All(p => p.HasDefaultValue));
                    return ctor.Invoke(ctor.GetParameters().Select(p => p.DefaultValue).ToArray());
                });

            foreach (var group in models.GroupBy(PolicyKey))
            {
                var model = group.First();
                var acceptsAll = (bool)model.GetType().GetProperty("AcceptsAllUnimodModifications")!.GetValue(model)!;
                var allowedIds = (IReadOnlySet<int>)model.GetType().GetProperty("AllowedUnimodIds")!.GetValue(model)!;
                if (!acceptsAll && allowedIds.Count == 0)
                    continue;
                var allowed = allowedIds.Where(unimodMasses.ContainsKey).Select(id => unimodMasses[id]).ToList();
                Assert.That(acceptsAll || allowed.Count > 0, Is.True, $"no UNIMOD masses found for {model.GetType().Name}'s allowed ids");
                yield return (model, acceptsAll ? _ => true : mass => allowed.Any(a => Math.Abs(a - mass) <= 2 * MassTolerance));
            }
        }

        private static string PolicyKey(object model)
        {
            var type = model.GetType();
            var acceptsAll = (bool)type.GetProperty("AcceptsAllUnimodModifications")!.GetValue(model)!;
            var allowed = (IReadOnlySet<int>)type.GetProperty("AllowedUnimodIds")!.GetValue(model)!;
            var required = (IReadOnlySet<int>?)type.GetProperty("RequiredNTerminalUnimodIds")!.GetValue(model);
            return $"{(acceptsAll ? "all" : string.Join(",", allowed.Order()))}|{(required == null ? "-" : string.Join(",", required.Order()))}";
        }

        // The checks are independent and the lookups' caches are concurrent, so they run in parallel to keep the
        // fixture fast.
        private static void AssertRoundTrips(Func<PeptideWithSetModifications, string> write, ISequenceParser? parser, int minimumAccepted,
            bool asWritten = false, bool validatedAsWritten = false)
        {
            var failures = new ConcurrentBag<string>();
            int accepted = 0;
            var work = ModelsByPolicy()
                .SelectMany(policy => Digests.Value.Where(d => policy.CouldAccept(d.ModificationMass)).Select(d => (policy.Model, d.Peptide)))
                .ToList();

            Parallel.ForEach(work, item =>
            {
                var (model, peptide) = item;
                var name = model.GetType().Name;
                var sequence = write(peptide);
                var validated = Clean(model, sequence, parser, out var payload);
                if (validated == null || payload == null)
                    return;
                Interlocked.Increment(ref accepted);
                if ((asWritten || validatedAsWritten) && validated != sequence)
                    failures.Add($"{name}: {sequence} validated as {validated}");

                var again = Clean(model, validated, parser, out var payloadAgain);
                if (again == null || payloadAgain != payload)
                {
                    failures.Add($"{name}: {sequence} -> {validated} sends {payload}, re-cleaned sends {payloadAgain ?? "nothing"}");
                    return;
                }

                foreach (var (source, built) in new[] { ("validated", Build(model, validated, parser)), ("input", Build(model, sequence, parser)) })
                {
                    if (built.Peptide == null)
                        failures.Add($"{name}: {sequence} -> {validated} sends {payload}, but its {source} sequence builds no peptide: {built.Warning}");
                    else if (Math.Abs(built.Peptide.MonoisotopicMass - peptide.MonoisotopicMass) > MassTolerance)
                        failures.Add($"{name}: {sequence} ({source}) built with mass off by {built.Peptide.MonoisotopicMass - peptide.MonoisotopicMass:F4}");
                    else if (!built.Peptide.AllModsOneIsNterminus.Keys.Order().SequenceEqual(peptide.AllModsOneIsNterminus.Keys.Order()))
                        failures.Add($"{name}: {sequence} ({source}) built with modifications at {string.Join(",", built.Peptide.AllModsOneIsNterminus.Keys.Order())} " +
                            $"instead of {string.Join(",", peptide.AllModsOneIsNterminus.Keys.Order())}");
                    else if (asWritten && built.Peptide.FullSequence != sequence)
                        failures.Add($"{name}: {sequence} ({source}) labeled {built.Peptide.FullSequence}");
                }
            });

            TestContext.Out.WriteLine($"{Digests.Value.Count} digested peptides, {work.Count} checks, {accepted} accepted, {failures.Count} failures");
            foreach (var line in failures.Order().Take(500))
                TestContext.Out.WriteLine(line);
            Assert.That(accepted, Is.GreaterThanOrEqualTo(minimumAccepted), "the corpus should exercise the models");
            Assert.That(failures, Is.Empty);
        }

        private static (PeptideWithSetModifications? Peptide, string? Warning) Build(object model, string sequence, ISequenceParser? parser)
        {
            var args = new object?[] { sequence, parser, null };
            var peptide = (PeptideWithSetModifications?)typeof(FragmentIntensityModel)
                .GetMethod("TryBuildPeptide", BindingFlags.NonPublic | BindingFlags.Instance)!.Invoke(model, args);
            return (peptide, ((WarningException?)args[2])?.Message);
        }

        private static string? Clean(object model, string sequence, ISequenceParser? parser, out string? payload)
        {
            var type = model.GetType();
            var tryClean = type.GetMethod("TryCleanSequence", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var args = new object?[] { sequence, parser, null, null };
            var validated = (string?)tryClean.Invoke(model, args);
            var koinaSequence = (CanonicalSequence?)args[2];
            payload = null;
            if (validated != null && koinaSequence != null)
            {
                var serialize = type.GetMethod("SerializeKoinaSequence", BindingFlags.NonPublic | BindingFlags.Instance)!;
                payload = (string?)serialize.Invoke(model, new object?[] { koinaSequence.Value, null });
            }
            return validated;
        }
    }
}
