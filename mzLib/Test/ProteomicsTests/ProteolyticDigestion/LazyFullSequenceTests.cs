using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;
using Chemistry;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Omics.SequenceConversion;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// A form made by digestion builds its full sequence when it is first read, not when it is made. These tests hold
    /// that full sequence to the one the previous code built in the constructor - kept verbatim below as the expression
    /// it evaluated, and also produced by the public constructor, which still builds it at once - and cover the three
    /// ways deferring it could change what a caller sees: binary serialization of the peptide index, ConvertModifications
    /// changing the modifications in place, and threads reading it for the first time together.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class LazyFullSequenceTests
    {
        private static readonly FieldInfo FullSequenceField = typeof(PeptideWithSetModifications)
            .GetField("<FullSequence>k__BackingField", BindingFlags.Instance | BindingFlags.NonPublic);

        private static string Stored(PeptideWithSetModifications form) => (string)FullSequenceField.GetValue(form);

        // the previous constructor's line: FullSequence = this.DetermineFullSequence();
        private static string FullSequenceBefore(PeptideWithSetModifications form) =>
            IBioPolymerWithSetMods.DetermineFullSequence(form.BaseSequence, form.AllModsOneIsNterminus);

        [Test]
        public static void DigestedFormsBuildTheFullSequenceThePreviousCodeBuilt()
        {
            var random = new Random(4);
            long forms = 0, modified = 0;
            for (int trial = 0; trial < 300; trial++)
            {
                DigestionParams digestionParams = RandomDigestionInputs.RandomDigestionParams(random);
                (List<Modification> fixedMods, List<Modification> variableMods) = RandomDigestionInputs.RandomSearchMods(random);
                Protein protein = RandomDigestionInputs.RandomProtein(random, trial, 5, 300, isDecoy: random.Next(4) == 0);

                foreach (PeptideWithSetModifications form in protein.Digest(digestionParams, fixedMods, variableMods))
                {
                    Assert.That(Stored(form), Is.Null, "a digested form must not build its full sequence before it is read");

                    var reference = new PeptideWithSetModifications(protein, digestionParams, form.OneBasedStartResidue, form.OneBasedEndResidue,
                        form.CleavageSpecificityForFdrCategory, form.PeptideDescription, form.MissedCleavages,
                        new Dictionary<int, Modification>(form.AllModsOneIsNterminus), form.NumFixedMods);
                    Assert.That(Stored(reference), Is.Not.Null, "the public constructor builds it at once, as before");

                    string expected = FullSequenceBefore(form);
                    Assert.That(form.FullSequence, Is.EqualTo(expected));
                    Assert.That(form.FullSequence, Is.EqualTo(reference.FullSequence));
                    Assert.That(Stored(form), Is.SameAs(form.FullSequence), "built once, then kept");
                    Assert.That(form.GetHashCode(), Is.EqualTo(reference.GetHashCode()));
                    Assert.That(form, Is.EqualTo(reference));
                    forms++;
                    if (form.AllModsOneIsNterminus.Count > 0)
                    {
                        modified++;
                    }
                }
            }
            Assert.That(forms, Is.GreaterThan(10000));
            Assert.That(modified, Is.GreaterThan(5000));
        }

        /// <summary>
        /// The peptide index cache is written with NetSerializer, field by field. Its layout must not change, or an index
        /// cached by the previous version would be misread: these are the serialized fields of mzLib 1.0.594.
        /// </summary>
        [Test]
        public static void SerializedFieldsAreThoseOfThePreviousVersion()
        {
            var expected = new[]
            {
                "DigestionProduct\t<CleavageSpecificityForFdrCategory>k__BackingField\tOmics.Digestion.CleavageSpecificity",
                "DigestionProduct\t<Description>k__BackingField\tSystem.String",
                "DigestionProduct\t<MissedCleavages>k__BackingField\tSystem.Int32",
                "DigestionProduct\t<OneBasedEndResidue>k__BackingField\tSystem.Int32",
                "DigestionProduct\t<OneBasedStartResidue>k__BackingField\tSystem.Int32",
                "DigestionProduct\t_baseSequence\tSystem.String",
                "PeptideWithSetModifications\t<FullSequence>k__BackingField\tSystem.String",
                "PeptideWithSetModifications\t<NumFixedMods>k__BackingField\tSystem.Int32",
                "PeptideWithSetModifications\t<PairedTargetDecoySequence>k__BackingField\tSystem.String",
                "PeptideWithSetModifications\tProteinAccession\tSystem.String",
                "PeptideWithSetModifications\t_fullSequenceWithMassShifts\tSystem.String",
            };

            var chain = new List<Type>();
            for (Type type = typeof(PeptideWithSetModifications); type != null && type != typeof(object); type = type.BaseType)
            {
                chain.Insert(0, type);
            }
            var actual = chain.SelectMany(type => type
                    .GetFields(BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.DeclaredOnly)
                    .Where(f => (f.Attributes & FieldAttributes.NotSerialized) == 0)
                    .OrderBy(f => f.Name, StringComparer.Ordinal)
                    .Select(f => type.Name + "\t" + f.Name + "\t" + f.FieldType.FullName))
                .ToArray();

            Assert.That(actual, Is.EqualTo(expected));
        }

        /// <summary>
        /// Digested forms written once their full sequence has been read - as MetaMorpheus writes its peptide index - come
        /// back as they went in. One written without reading it carries no full sequence, which is why the writer reads it.
        /// </summary>
        [Test]
        public static void DigestedFormsRoundTripThroughBinarySerializationOnceRead()
        {
            var random = new Random(9);
            var digestionParams = new DigestionParams("trypsin", 2, 5);
            var proteins = Enumerable.Range(0, 20).Select(i => RandomDigestionInputs.RandomProtein(random, i, 30, 200)).ToList();
            var fixedMods = new List<Modification> { RandomDigestionInputs.CarbamidomethylC };
            var variableMods = new List<Modification> { RandomDigestionInputs.OxidationM };
            List<PeptideWithSetModifications> written = proteins.SelectMany(p => p.Digest(digestionParams, fixedMods, variableMods)).ToList();
            Assert.That(written.Count(p => p.AllModsOneIsNterminus.Count > 1), Is.GreaterThan(50));

            PeptideWithSetModifications neverRead = written[0];
            var unread = (PeptideWithSetModifications)RoundTrip(new List<PeptideWithSetModifications> { neverRead }, "unread")[0];
            Assert.That(Stored(unread), Is.Null, "nothing to write until the full sequence is read");

            foreach (PeptideWithSetModifications form in written)
            {
                _ = form.FullSequence;
            }
            List<PeptideWithSetModifications> read = RoundTrip(written, "read");

            var idToMod = RandomDigestionInputs.Localized.Concat(fixedMods).Concat(variableMods).Distinct().ToDictionary(m => m.IdWithMotif);
            var accessionToProtein = proteins.ToDictionary(p => p.Accession);
            Assert.That(read.Count, Is.EqualTo(written.Count));
            for (int i = 0; i < written.Count; i++)
            {
                read[i].SetNonSerializedPeptideInfo(idToMod, accessionToProtein, digestionParams);
                Assert.That(read[i].FullSequence, Is.EqualTo(written[i].FullSequence));
                Assert.That(read[i], Is.EqualTo(written[i]));
                Assert.That(read[i].MonoisotopicMass, Is.EqualTo(written[i].MonoisotopicMass).Within(1e-9));
                Assert.That(read[i].AllModsOneIsNterminus.OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value)),
                    Is.EqualTo(written[i].AllModsOneIsNterminus.OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value))));
            }
        }

        private static List<PeptideWithSetModifications> RoundTrip(List<PeptideWithSetModifications> peptides, string name)
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "LazyFullSequenceTests");
            Directory.CreateDirectory(dir);
            string path = Path.Combine(dir, name + ".ind");
            var serializer = new NetSerializer.Serializer(new List<Type> { typeof(List<PeptideWithSetModifications>), typeof(PeptideWithSetModifications) });
            using (var file = File.Create(path))
            {
                serializer.Serialize(file, peptides);
            }
            using (var file = File.OpenRead(path))
            {
                return (List<PeptideWithSetModifications>)serializer.Deserialize(file);
            }
        }

        /// <summary>
        /// ConvertModifications changes a form's modifications in place. The full sequence was built before that, in the
        /// constructor, and must still be: here the lookup resolves nothing, so every modification is removed.
        /// </summary>
        [Test]
        public static void ConvertingModificationsLeavesTheFullSequenceTheFormWasMadeWith()
        {
            var random = new Random(5);
            int converted = 0;
            for (int trial = 0; trial < 50; trial++)
            {
                Protein protein = RandomDigestionInputs.RandomProtein(random, trial, 20, 120);
                var digestionParams = new DigestionParams("trypsin", 1, 5);
                foreach (PeptideWithSetModifications form in protein.Digest(digestionParams, new List<Modification> { RandomDigestionInputs.CarbamidomethylC },
                             new List<Modification> { RandomDigestionInputs.OxidationM }).Where(f => f.AllModsOneIsNterminus.Count > 0))
                {
                    var reference = new PeptideWithSetModifications(protein, digestionParams, form.OneBasedStartResidue, form.OneBasedEndResidue,
                        form.CleavageSpecificityForFdrCategory, form.PeptideDescription, form.MissedCleavages,
                        new Dictionary<int, Modification>(form.AllModsOneIsNterminus), form.NumFixedMods);
                    Assert.That(Stored(form), Is.Null);

                    form.ConvertModifications(new ResolvesNothing());
                    reference.ConvertModifications(new ResolvesNothing());

                    Assert.That(form.AllModsOneIsNterminus, Is.Empty, "the conversion must have changed the modifications");
                    Assert.That(form.FullSequence, Is.EqualTo(reference.FullSequence));
                    Assert.That(form.FullSequence, Does.Contain("["));
                    converted++;
                }
            }
            Assert.That(converted, Is.GreaterThan(100));
        }

        private class ResolvesNothing : IModificationLookup
        {
            public string Name => "resolves nothing";
            public CanonicalModification? TryResolve(CanonicalModification mod) => null;
            public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null,
                ChemicalFormula chemicalFormula = null, ModificationPositionType? positionType = null) => null;
        }

        [Test]
        public static void ThreadsReadingTheFullSequenceFirstTogetherAllGetIt()
        {
            var random = new Random(6);
            Protein protein = RandomDigestionInputs.RandomProtein(random, 0, 300, 300);
            var digestionParams = new DigestionParams("top-down", maxModsForPeptides: 3);
            List<PeptideWithSetModifications> forms = protein.Digest(digestionParams, new List<Modification>(),
                new List<Modification> { RandomDigestionInputs.OxidationM }).ToList();
            Assert.That(forms.Count, Is.GreaterThan(100));

            // dedicated threads, released together, so the first reads really do overlap
            const int threads = 8;
            var seen = new string[threads][];
            using var barrier = new Barrier(threads);
            var readers = Enumerable.Range(0, threads).Select(t => new Thread(() =>
            {
                barrier.SignalAndWait();
                seen[t] = forms.Select(f => f.FullSequence).ToArray();
            })).ToList();
            readers.ForEach(r => r.Start());
            readers.ForEach(r => r.Join());

            for (int i = 0; i < forms.Count; i++)
            {
                string expected = FullSequenceBefore(forms[i]);
                for (int t = 0; t < threads; t++)
                {
                    Assert.That(seen[t][i], Is.EqualTo(expected));
                }
                Assert.That(forms[i].FullSequence, Is.EqualTo(expected));
            }
        }
    }
}
