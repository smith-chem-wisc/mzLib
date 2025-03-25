using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Diagnostics.CodeAnalysis;
using Transcriptomics;
using Transcriptomics.Digestion;
using System.IO;
using NetSerializer;
using MzLibUtil;
using NUnit.Framework.Legacy;

namespace Test.FileReadingTests;

[ExcludeFromCodeCoverage]
public class TestPeptideSerializer
{
    private static void RefreshDirectory(string path)
    {
        if (Directory.Exists(path))
        {
            Directory.Delete(path, true);
        }
        Directory.CreateDirectory(path!);
    }

    // Peptide tests that were moved and no changes made to the tests
    [Test]
    public static void TestSerializationPeptideFromString()
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is unmodified and generated from reading in a string
        string sequence = "PEPTIDE";
        PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
        PeptideWithSetModifications deserializedPeptide = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromString");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");

        var messageTypes = typeof(PeptideWithSetModifications);
        var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }

        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer>(), null);

        // not asserting any rna properties - since the oligo was created from a sequence string it didn't have a rna to begin with

        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

        var products = new List<Product>();

        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }

    [Test]
    public static void TestSerializationPeptideFromProtein()
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is unmodified and generated from digesting a rna
        Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein");

        PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
        PeptideWithSetModifications deserializedPeptide = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");

        var messageTypes = typeof(PeptideWithSetModifications);
        var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }

        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

        Assert.That(peptide.DigestionParams.Equals(deserializedPeptide.DigestionParams));
        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

        var products = new List<Product>();

        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }

    [Test]
    public static void TestSerializationPeptideFromProteinWithMod()
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is modified with a phosphorylation

        ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

        Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
            };

        Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

        Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

        Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein", oneBasedModifications: mods);

        PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.AllModsOneIsNterminus.Count == 1).First();
        PeptideWithSetModifications deserializedPeptide = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProteinWithMod");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");

        var messageTypes = typeof(PeptideWithSetModifications);
        var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }

        Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };

        deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

        var products = new List<Product>();

        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }

    // Peptide Tests with legacy and serializable sequence interface
    [Test]
    [TestCase(1, TestName ="All Manually")]
    [TestCase(2, TestName ="Serializer Manually - Types from Interface")]
    [TestCase(3, TestName ="All From Interface")]
    public static void TestSerializationPeptideFromString_Interface(int serializerConstructionMethod)
    {
        string sequence = "PEPTIDE";
        PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
        PeptideWithSetModifications deserializedPeptide = null;
        
        // Load serializer
        Serializer ser;
        Type[] messagingTypes;
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromString_new");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        switch (serializerConstructionMethod)
        {
            case 1:
                messagingTypes = new Type[] { typeof(PeptideWithSetModifications) };
                ser = new Serializer(messagingTypes);
                break;

            case 2:
                messagingTypes = peptide.GetTypesToSerialize();
                ser = new Serializer(messagingTypes);
                break;

            case 3:
                ser = peptide.GetSequenceSerializer();
                break;

            default:
                throw new MzLibException("Test Case Not Implemented");
        }

        // Serialize
        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        // Deserialize
        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }
        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer>(), null);

        // not asserting any rna properties - since the oligo was created from a sequence string it didn't have a rna to begin with
        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

        var products = new List<Product>();
        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
    }

    [Test]
    [TestCase(1, TestName = "All Manually")]
    [TestCase(2, TestName = "Serializer Manually - Types from Interface")]
    [TestCase(3, TestName = "All From Interface")]
    public static void TestSerializationPeptideFromProtein_Interface(int serializerConstructionMethod)
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is unmodified and generated from digesting a rna
        Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein");

        PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
        PeptideWithSetModifications deserializedPeptide = null;

        // Load serializer
        Serializer ser;
        Type[] messagingTypes;
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein_new");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        switch (serializerConstructionMethod)
        {
            case 1:
                messagingTypes = new Type[] { typeof(PeptideWithSetModifications) };
                ser = new Serializer(messagingTypes);
                break;

            case 2:
                messagingTypes = peptide.GetTypesToSerialize();
                ser = new Serializer(messagingTypes);
                break;

            case 3:
                ser = peptide.GetSequenceSerializer();
                break;

            default:
                throw new MzLibException("Test Case Not Implemented");
        }

        // Serialize
        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        // Deserialize
        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }
        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

        Assert.That(peptide.DigestionParams.Equals(deserializedPeptide.DigestionParams));
        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

        var products = new List<Product>();
        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
    }

    [Test]
    [TestCase(1, TestName = "All Manually")]
    [TestCase(2, TestName = "Serializer Manually - Types from Interface")]
    [TestCase(3, TestName = "All From Interface")]
    public static void TestSerializationPeptideFromProteinWithMod_Interface(int serializerConstructionMethod)
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is modified with a phosphorylation
        ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

        Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
        {
            { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
            { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
        };

        Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

        Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

        Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein", oneBasedModifications: mods);

        PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.AllModsOneIsNterminus.Count == 1).First();
        PeptideWithSetModifications deserializedPeptide = null;

        // Load serializer
        Serializer ser;
        Type[] messagingTypes;
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein_new");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        switch (serializerConstructionMethod)
        {
            case 1:
                messagingTypes = new Type[] { typeof(PeptideWithSetModifications) };
                ser = new Serializer(messagingTypes);
                break;

            case 2:
                messagingTypes = peptide.GetTypesToSerialize();
                ser = new Serializer(messagingTypes);
                break;

            case 3:
                ser = peptide.GetSequenceSerializer();
                break;

            default:
                throw new MzLibException("Test Case Not Implemented");
        }

        // Serialize
        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        // Deserialize
        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
        }
        Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };
        deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);
        Assert.That(deserializedPeptide.AllModsOneIsNterminus.Count, Is.EqualTo(peptide.AllModsOneIsNterminus.Count));
        CollectionAssert.AreEqual(deserializedPeptide.AllModsOneIsNterminus, peptide.AllModsOneIsNterminus);

        var products = new List<Product>();
        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }

    [Test]
    [TestCase(1, TestName = "All Manually")]
    [TestCase(2, TestName = "Serializer Manually - Types from Interface")]
    [TestCase(3, TestName = "All From Interface")]
    public static void TestSerializationPeptideFromProteinWithMod_InterfaceAsIBioPolymerWithSetMods(int serializerConstructionMethod)
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is modified with a phosphorylation
        ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

        Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
        {
            { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
            { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
        };

        Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

        Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

        Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein", oneBasedModifications: mods);

        IBioPolymerWithSetMods peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First(v => v.AllModsOneIsNterminus.Count == 1);
        IBioPolymerWithSetMods deserializedPeptide = null;

        // Load serializer
        Serializer ser;
        Type[] messagingTypes;
        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein_new");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        switch (serializerConstructionMethod)
        {
            case 1:
                messagingTypes = new Type[] { typeof(PeptideWithSetModifications) };
                ser = new Serializer(messagingTypes);
                break;

            case 2:
                messagingTypes = peptide.GetTypesToSerialize();
                ser = new Serializer(messagingTypes);
                break;

            case 3:
                ser = peptide.GetSequenceSerializer();
                break;

            default:
                throw new MzLibException("Test Case Not Implemented");
        }

        // Serialize
        using (var file = File.Create(path))
        {
            ser.Serialize(file, peptide);
        }

        // Deserialize
        using (var file = File.OpenRead(path))
        {
            deserializedPeptide = (IBioPolymerWithSetMods)ser.Deserialize(file);
        }
        Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };
        deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

        Assert.That(peptide.Equals(deserializedPeptide));
        Assert.That(deserializedPeptide.Parent.Name == peptide.Parent.Name);
        Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
        Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);
        Assert.That(deserializedPeptide.AllModsOneIsNterminus.Count, Is.EqualTo(peptide.AllModsOneIsNterminus.Count));
        CollectionAssert.AreEqual(deserializedPeptide.AllModsOneIsNterminus, peptide.AllModsOneIsNterminus);

        var products = new List<Product>();
        deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }


    // Oligo Tests that are a copy of the above oligo tests
    [Test]
    public static void TestSerializationOligoFromString()
    {
        // purpose of this test is to serialize/deserialize a OligoWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is unmodified and generated from reading in a string
        string sequence = "GUACUGAGUCUACUAGAUCA";
        OligoWithSetMods oligo = new OligoWithSetMods(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
        OligoWithSetMods deserializeOligo = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerialization_OligoWithSetMods_FromString");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myOligoIndex.ind");

        var ser = new NetSerializer.Serializer(new List<Type> { typeof(OligoWithSetMods), typeof(ChemicalFormula) });
        var typeMap = ser.GetTypeMap();

        using (var file = File.Create(path))
        {
            ser.Serialize(file, oligo);
        }

        using (var file = File.OpenRead(path))
        {
            deserializeOligo = (OligoWithSetMods)ser.Deserialize(file);
        }

        deserializeOligo.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer>(), null);

        // not asserting any rna properties - since the oligo was created from a sequence string it didn't have a rna to begin with

        Assert.That(oligo.Equals(deserializeOligo), Is.True);
        Assert.That(deserializeOligo.MonoisotopicMass, Is.EqualTo(oligo.MonoisotopicMass));
        Assert.That(deserializeOligo.SequenceWithChemicalFormulas, Is.EqualTo(oligo.SequenceWithChemicalFormulas));
        Assert.That(deserializeOligo.FivePrimeTerminus, Is.EqualTo(oligo.FivePrimeTerminus));
        Assert.That(deserializeOligo.ThreePrimeTerminus, Is.EqualTo(oligo.ThreePrimeTerminus));

        var products = new List<Product>();

        deserializeOligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        oligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }

    [Test]
    public static void TestSerializationOligoFromRna()
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is unmodified and generated from digesting a rna
        RNA rna = new("GUACUGAGUCUACUAGAUCA", "name", "Accession1", "", "");

        OligoWithSetMods oligo = rna.Digest(new RnaDigestionParams("RNase T1"), new List<Modification>(), new List<Modification>()).First();
        OligoWithSetMods deserializeOligo = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerialization_OligoWithSetMods_FromRNA");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        var ser = new Serializer(new List<Type> { typeof(OligoWithSetMods), typeof(ChemicalFormula) });

        using (var file = File.Create(path))
        {
            ser.Serialize(file, oligo);
        }

        using (var file = File.OpenRead(path))
        {
            deserializeOligo = (OligoWithSetMods)ser.Deserialize(file);
        }
        deserializeOligo.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer> { { rna.Accession, rna } }, oligo.DigestionParams);

        Assert.That(oligo.DigestionParams.Equals(deserializeOligo.DigestionParams));
        Assert.That(deserializeOligo.Parent.Name == oligo.Parent.Name);
        Assert.That(oligo.Equals(deserializeOligo), Is.True);
        Assert.That(deserializeOligo.MonoisotopicMass, Is.EqualTo(oligo.MonoisotopicMass));
        Assert.That(deserializeOligo.SequenceWithChemicalFormulas, Is.EqualTo(oligo.SequenceWithChemicalFormulas));
        Assert.That(deserializeOligo.FivePrimeTerminus, Is.EqualTo(oligo.FivePrimeTerminus));
        Assert.That(deserializeOligo.ThreePrimeTerminus, Is.EqualTo(oligo.ThreePrimeTerminus));

        var products = new List<Product>();

        deserializeOligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        oligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }



    [Test]
    public static void TestSerializationOligoFromRnaWithMod()
    {
        // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized oligo
        // has the same properties as before it was serialized. This oligo is modified with a phosphorylation

        ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

        Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
        {
            { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
            { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
        };

        Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

        Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

        RNA rna = new("GUACUGAGUCUACUAGAUCA", "name", "Accession1", "", "");

        OligoWithSetMods oligo = rna.Digest(new RnaDigestionParams("RNase T1"), new List<Modification>(), new List<Modification>()).First();
        OligoWithSetMods deserializeOligo = null;

        string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerialization_OligoWithSetMods_FromRNA_WithMod");
        RefreshDirectory(dir);
        string path = Path.Combine(dir, "myPeptideIndex.ind");
        var ser = new Serializer(new List<Type> { typeof(OligoWithSetMods), typeof(ChemicalFormula) });

        using (var file = File.Create(path))
        {
            ser.Serialize(file, oligo);
        }

        using (var file = File.OpenRead(path))
        {
            deserializeOligo = (OligoWithSetMods)ser.Deserialize(file);
        }

        Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };

        deserializeOligo.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { rna.Accession, rna } }, oligo.DigestionParams);

        Assert.That(oligo.Equals(deserializeOligo));
        Assert.That(deserializeOligo.Parent.Name == oligo.Parent.Name);
        Assert.That(deserializeOligo.MonoisotopicMass, Is.EqualTo(oligo.MonoisotopicMass));
        Assert.That(deserializeOligo.SequenceWithChemicalFormulas, Is.EqualTo(oligo.SequenceWithChemicalFormulas));
        Assert.That(deserializeOligo.FivePrimeTerminus, Is.EqualTo(oligo.FivePrimeTerminus));
        Assert.That(deserializeOligo.ThreePrimeTerminus, Is.EqualTo(oligo.ThreePrimeTerminus));

        var products = new List<Product>();

        deserializeOligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

        oligo.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

        Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        Directory.Delete(dir, true);
    }
}