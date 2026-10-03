using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// ValidatedFullSequence and the Koina payload must name the same modifications. Over every UNIMOD modification
    /// written as a ProForma mass and as a ProForma name, and every mzLib modification name, at residues and termini, on
    /// one model per distinct modification policy: re-cleaning an accepted ValidatedFullSequence with the same parser
    /// gives the same payload, and the peptide built from its modifications has the mass, and the modifications at the
    /// same places, as the peptide built from the payload's. mzLib names mzLib knows label that peptide as written.
    /// </summary>
    [TestFixture]
    public class KoinaValidatedSequenceRoundTripTests
    {
        private const string Residues = "ACDEFGHIKLMNPQRSTVWY";
        private const double MassTolerance = 0.01;

        [Test]
        public void ValidatedFullSequence_RoundTripsToThePayloadAndItsMass_ProFormaMasses()
        {
            var corpus = Corpus(Mods.UnimodModifications, m => m.MonoisotopicMass!.Value.ToString("+0.0000;-0.0000", CultureInfo.InvariantCulture), "-");

            // Without a residue and position check, the models' lookup matches a C-terminal -18.0106 on Q or N to
            // N-terminal pyro-Glu (UNIMOD:27); its id, which is all ValidatedFullSequence carries, then finds nothing
            // at the C-terminus. The lookup's mismatch, not a disagreement, so these and only these are allowed.
            AssertRoundTrips(corpus, ProFormaSequenceParser.Instance, minimumAccepted: 1000, expectedMisfits: new HashSet<string>
            {
                "Prosit2024IntensityPTMsGl: GLSTDEFN-[-18.0106]", "Prosit2024IntensityPTMsGl: GLSTDEFQ-[-18.0106]",
                "UniSpec: GLSTDEFN-[-18.0106]", "UniSpec: GLSTDEFQ-[-18.0106]"
            });
        }

        [Test]
        public void ValidatedFullSequence_RoundTripsToThePayloadAndItsMass_ProFormaNames()
        {
            var corpus = Corpus(Mods.UnimodModifications.Where(m => !string.IsNullOrWhiteSpace(m.OriginalId) && !m.OriginalId.Contains('[') && !m.OriginalId.Contains(']')),
                m => m.OriginalId, "-");

            AssertRoundTrips(corpus, ProFormaSequenceParser.Instance, minimumAccepted: 1000);
        }

        [Test]
        public void ValidatedFullSequence_RoundTripsToThePayloadAndItsMass_MzLibNames()
        {
            var corpus = Corpus(Mods.AllKnownProteinModsDictionary.Values.Where(m => !string.IsNullOrEmpty(m.ModificationType)),
                m => $"{m.ModificationType}:{m.IdWithMotif}", "");

            AssertRoundTrips(corpus, null, minimumAccepted: 1000);
        }

        [Test]
        public void ValidatedFullSequence_ProFormaName_KeepsItsText()
        {
            var model = new Ms2PipHCD2021();

            var validated = Clean(model, "GLSK[Acetyl]DEFK", ProFormaSequenceParser.Instance, out var payload, out _);

            Assert.That(validated, Is.EqualTo("GLSK[Acetyl]DEFK"));
            Assert.That(payload, Does.Contain("UNIMOD:1]"));
        }

        [TestCase("GLST[-1.9793]DEFK", 1210)]
        [TestCase("GLSS[+79.9568]DEFK", 21)]
        public void ValidatedFullSequence_ProFormaMassOnlyModification_NamesWhatKoinaIsSent(string sequence, int unimodId)
        {
            var model = new Ms2PipHCD2021();

            var validated = Clean(model, sequence, ProFormaSequenceParser.Instance, out var payload, out _);

            Assert.That(payload, Does.Contain($"UNIMOD:{unimodId}]"));
            Assert.That(validated, Is.EqualTo(payload!.Replace("unimod", "UNIMOD")));
            Assert.That(Clean(model, validated!, ProFormaSequenceParser.Instance, out var again, out _), Is.Not.Null);
            Assert.That(again, Is.EqualTo(payload));
        }

        // One peptide per modification and place: on its residue, at the N-terminus or at the C-terminus.
        private static List<(string Sequence, double Mass)> Corpus(IEnumerable<Modification> mods, Func<Modification, string> text, string terminalSeparator) =>
            mods.Where(m => m.MonoisotopicMass.HasValue && m.Target != null && m.Target.ToString().Length == 1
                    && (Residues.Contains(m.Target.ToString()) || m.Target.ToString() == "X"))
                .Select(m =>
                {
                    var target = m.Target.ToString();
                    var mod = $"[{text(m)}]";
                    string? sequence = m.LocationRestriction switch
                    {
                        "Anywhere." when target != "X" => $"GLS{target}{mod}DEFK",
                        "N-terminal." or "Peptide N-terminal." => $"{mod}{terminalSeparator}{(target == "X" ? "G" : target)}LSTDEFK",
                        "C-terminal." or "Peptide C-terminal." => $"GLSTDEF{(target == "X" ? "K" : target)}-{mod}",
                        _ => null
                    };
                    return (Sequence: sequence, Mass: m.MonoisotopicMass!.Value);
                })
                .Where(entry => entry.Sequence != null)
                .Select(entry => (entry.Sequence!, entry.Mass))
                .DistinctBy(entry => entry.Item1)
                .ToList();

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
                    continue; // accepts no modification at all
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
        private static void AssertRoundTrips(List<(string Sequence, double Mass)> corpus, ISequenceParser? parser, int minimumAccepted,
            IReadOnlySet<string>? expectedMisfits = null)
        {
            var failures = new ConcurrentBag<string>();
            var misfits = new ConcurrentBag<string>();
            int accepted = 0;
            var work = ModelsByPolicy()
                .SelectMany(policy => corpus.Where(entry => policy.CouldAccept(entry.Mass)).Select(entry => (policy.Model, entry.Sequence)))
                .ToList();

            Parallel.ForEach(work, item =>
            {
                var (model, sequence) = item;
                var name = model.GetType().Name;
                var validated = Clean(model, sequence, parser, out var payload, out var koinaSequence);
                if (validated == null || payload == null)
                    return;
                Interlocked.Increment(ref accepted);

                var again = Clean(model, validated, parser, out var payloadAgain, out _);
                if (again == null || payloadAgain != payload)
                {
                    failures.Add($"{name}: {sequence} -> {validated} sends {payload}, re-cleaned sends {payloadAgain ?? "nothing"}");
                    return;
                }

                var fromPayload = Build(koinaSequence!.Value.BaseSequence, koinaSequence.Value.Modifications.Select(m => (m, m.MzLibModification!)));
                var parsed = (parser ?? MzLibSequenceParser.Instance).Parse(validated, null, SequenceConversionHandlingMode.ThrowException)!.Value;
                var resolved = parsed.Modifications.Select(m => (m, Resolve(model, m))).ToList();
                if (resolved.Any(m => m.Item2 == null))
                {
                    if (expectedMisfits?.Contains($"{name}: {sequence}") == true)
                        misfits.Add($"{name}: {sequence}");
                    else
                        failures.Add($"{name}: {sequence} -> {validated} sends {payload}, but a modification of it resolves to nothing");
                    return;
                }
                var fromValidated = Build(parsed.BaseSequence, resolved.Select(m => (m.Item1, m.Item2!)));

                if (Math.Abs(fromValidated.MonoisotopicMass - fromPayload.MonoisotopicMass) > MassTolerance)
                    failures.Add($"{name}: {sequence} -> {validated} sends {payload}, mass off by {fromValidated.MonoisotopicMass - fromPayload.MonoisotopicMass:F4}");
                else if (parser == null && fromValidated.FullSequence != validated)
                    failures.Add($"{name}: {sequence} -> {validated} is labeled {fromValidated.FullSequence}, not as written");
                else if (!fromValidated.AllModsOneIsNterminus.Keys.Order().SequenceEqual(fromPayload.AllModsOneIsNterminus.Keys.Order()))
                    failures.Add($"{name}: {sequence} -> {validated} sends {payload}, modifications at {string.Join(",", fromValidated.AllModsOneIsNterminus.Keys.Order())} " +
                        $"instead of {string.Join(",", fromPayload.AllModsOneIsNterminus.Keys.Order())}");
            });

            TestContext.Out.WriteLine($"{corpus.Count} inputs, {work.Count} checks, {accepted} accepted, {failures.Count} failures, " +
                $"{misfits.Count} expected lookup misfits");
            foreach (var line in failures.Order().Concat(misfits.Order()).Take(500))
                TestContext.Out.WriteLine(line);
            Assert.That(accepted, Is.GreaterThanOrEqualTo(minimumAccepted), "the corpus should exercise the models");
            Assert.That(failures, Is.Empty);
            Assert.That(misfits, Is.EquivalentTo(expectedMisfits ?? new HashSet<string>()), "a pinned lookup misfit no longer occurs");
        }

        private static Modification? Resolve(object model, CanonicalModification mod) =>
            (Modification?)typeof(FragmentIntensityModel).GetMethod("ResolveModificationObject", BindingFlags.NonPublic | BindingFlags.Instance)!
                .Invoke(model, new object?[] { mod });

        private static PeptideWithSetModifications Build(string baseSequence, IEnumerable<(CanonicalModification, Modification)> modifications) =>
            (PeptideWithSetModifications)typeof(FragmentIntensityModel)
                .GetMethod("BuildPeptide", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { baseSequence, modifications })!;

        private static string? Clean(object model, string sequence, ISequenceParser? parser, out string? payload, out CanonicalSequence? koinaSequence)
        {
            var type = model.GetType();
            var tryClean = type.GetMethod("TryCleanSequence", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var args = new object?[] { sequence, parser, null, null };
            var validated = (string?)tryClean.Invoke(model, args);
            koinaSequence = (CanonicalSequence?)args[2];
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
