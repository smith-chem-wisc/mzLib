using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Fragmentation;
using Transcriptomics.Digestion;

namespace Test.Transcriptomics;

public class BackboneFragmentationTestCase
{
    public string UnmodifiedSequence { get; }
    public string ModifiedSequence { get; }

    public Dictionary<ProductType, Product[]> UnmodifiedProducts { get; }
    public Dictionary<ProductType, Product[]> ModifiedProducts { get; }
    public Dictionary<ProductType, double[]> ExpectedDeltaMass { get; }

    public BackboneFragmentationTestCase(string originalSeq, string modifiedSeq, Dictionary<ProductType, double[]> expectedDelta)
    {
        UnmodifiedSequence = originalSeq;
        ModifiedSequence = modifiedSeq;
        ExpectedDeltaMass = expectedDelta;

        var unMod = new OligoWithSetMods(originalSeq);
        var mod = new OligoWithSetMods(modifiedSeq);

        var unModProd = new List<Product>();
        var modProd = new List<Product>();

        unMod.Fragment(MassSpectrometry.DissociationType.CID, FragmentationTerminus.Both, unModProd);
        mod.Fragment(MassSpectrometry.DissociationType.CID, FragmentationTerminus.Both, modProd);

        UnmodifiedProducts = unModProd.GroupBy(p => p.ProductType).ToDictionary(g => g.Key, g => g.ToArray());
        ModifiedProducts = modProd.GroupBy(p => p.ProductType).ToDictionary(g => g.Key, g => g.ToArray());
    }
}

[TestFixture]
public static class BackboneFragmentation
{
    public static IEnumerable<BackboneFragmentationTestCase> GetTestCases()
    {
        var modMass = ChemicalFormula.ParseFormula("SO-1").MonoisotopicMass;

        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUA[Common Variable:Phosphorothioate on X]CUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, 0, 0, modMass, modMass] },
                { ProductType.c,          [0, 0, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [0, 0, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, modMass, modMass, modMass] },
                { ProductType.y,          [0, 0, 0, modMass, modMass] },
                { ProductType.yWaterLoss, [0, 0, 0, modMass, modMass] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "G[Common Variable:Phosphorothioate on X]UACUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, modMass, modMass, modMass, modMass] },
                { ProductType.c,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, 0, 0, modMass] },
                { ProductType.y,          [0, 0, 0, 0, 0] },
                { ProductType.yWaterLoss, [0, 0, 0, 0, 0] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUACU[Common Variable:Phosphorothioate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, 0, 0] },
                { ProductType.aBaseLoss,  [0, 0, 0, 0, 0] },
                { ProductType.c,          [0, 0, 0, 0, modMass] },
                { ProductType.dWaterLoss, [0, 0, 0, 0, modMass] },
                { ProductType.w,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.y,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.yWaterLoss, [0, modMass, modMass, modMass, modMass] },
            }
        );

        modMass = ChemicalFormula.ParseFormula("B1 H3 O-1").MonoisotopicMass;

        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUA[Therapeutic:Boranophosphate on X]CUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, 0, 0, modMass, modMass] },
                { ProductType.c,          [0, 0, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [0, 0, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, modMass, modMass, modMass] },
                { ProductType.y,          [0, 0, 0, modMass, modMass] },
                { ProductType.yWaterLoss, [0, 0, 0, modMass, modMass] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "G[Therapeutic:Boranophosphate on X]UACUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, modMass, modMass, modMass, modMass] },
                { ProductType.c,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, 0, 0, modMass] },
                { ProductType.y,          [0, 0, 0, 0, 0] },
                { ProductType.yWaterLoss, [0, 0, 0, 0, 0] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUACU[Therapeutic:Boranophosphate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, 0, 0] },
                { ProductType.aBaseLoss,  [0, 0, 0, 0, 0] },
                { ProductType.c,          [0, 0, 0, 0, modMass] },
                { ProductType.dWaterLoss, [0, 0, 0, 0, modMass] },
                { ProductType.w,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.y,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.yWaterLoss, [0, modMass, modMass, modMass, modMass] },
            }
        );

        modMass = ChemicalFormula.ParseFormula("N1 H1 O-1").MonoisotopicMass;

        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUA[Therapeutic:N3'-Phosphoramidate on X]CUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, 0, 0, modMass, modMass] },
                { ProductType.c,          [0, 0, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [0, 0, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, modMass, modMass, modMass] },
                { ProductType.y,          [0, 0, 0, modMass, modMass] },
                { ProductType.yWaterLoss, [0, 0, 0, modMass, modMass] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "G[Therapeutic:N3'-Phosphoramidate on X]UACUG",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, modMass, modMass, modMass, modMass] },
                { ProductType.c,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, 0, 0, modMass] },
                { ProductType.y,          [0, 0, 0, 0, 0] },
                { ProductType.yWaterLoss, [0, 0, 0, 0, 0] },
            }
        );
        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "GUACU[Therapeutic:N3'-Phosphoramidate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, 0, 0] },
                { ProductType.aBaseLoss,  [0, 0, 0, 0, 0] },
                { ProductType.c,          [0, 0, 0, 0, modMass] },
                { ProductType.dWaterLoss, [0, 0, 0, 0, modMass] },
                { ProductType.w,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.y,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.yWaterLoss, [0, modMass, modMass, modMass, modMass] },
            }
        );


        modMass = ChemicalFormula.ParseFormula("SO-1").MonoisotopicMass;

        yield return new BackboneFragmentationTestCase(
            "GUACUG",
            "G[Common Variable:Phosphorothioate on X]UACU[Common Variable:Phosphorothioate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, modMass, modMass, modMass, modMass] },
                { ProductType.c,          [modMass, modMass, modMass, modMass, 2*modMass] },
                { ProductType.dWaterLoss, [modMass, modMass, modMass, modMass, 2*modMass] },
                { ProductType.w,          [modMass, modMass, modMass, modMass, 2*modMass] },
                { ProductType.y,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.yWaterLoss, [0, modMass, modMass, modMass, modMass] },
            }
        );

        yield return new BackboneFragmentationTestCase(
            "G[Common Variable:Phosphorothioate on X]UACUG",
            "G[Common Variable:Phosphorothioate on X]UACU[Common Variable:Phosphorothioate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, 0, 0, 0, 0] },
                { ProductType.aBaseLoss,  [0, 0, 0, 0, 0] },
                { ProductType.c,          [0, 0, 0, 0, modMass] },
                { ProductType.dWaterLoss, [0, 0, 0, 0, modMass] },
                { ProductType.w,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.y,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.yWaterLoss, [0, modMass, modMass, modMass, modMass] },
            }
        );

        yield return new BackboneFragmentationTestCase(
            "GUACU[Common Variable:Phosphorothioate on X]G",
            "G[Common Variable:Phosphorothioate on X]UACU[Common Variable:Phosphorothioate on X]G",
            new Dictionary<ProductType, double[]>
            {
                { ProductType.a,          [0, modMass, modMass, modMass, modMass] },
                { ProductType.aBaseLoss,  [0, modMass, modMass, modMass, modMass] },
                { ProductType.c,          [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.dWaterLoss, [modMass, modMass, modMass, modMass, modMass] },
                { ProductType.w,          [0, 0, 0, 0, modMass] },
                { ProductType.y,          [0, 0, 0, 0, 0] },
                { ProductType.yWaterLoss, [0, 0, 0, 0, 0] },
            }
        );
    }


    [Test]
    [TestCaseSource(nameof(GetTestCases))]
    public static void ProductsWithBackboneModificationHaveExpectedMassShifts(BackboneFragmentationTestCase testCase)
    {
        foreach (var productType in testCase.ExpectedDeltaMass.Keys)
        {
            var unmodifiedProducts = testCase.UnmodifiedProducts[productType];
            var modifiedProducts = testCase.ModifiedProducts[productType];
            var expectedDeltas = testCase.ExpectedDeltaMass[productType];
            for (int i = 0; i < expectedDeltas.Length; i++)
            {
                double actualDelta = modifiedProducts[i].MonoisotopicMass - unmodifiedProducts[i].MonoisotopicMass;
                Assert.That(actualDelta, Is.EqualTo(expectedDeltas[i]).Within(1e-6), $"Failed {testCase.ModifiedSequence} - {testCase.UnmodifiedSequence} for product type {productType}{i+1}");
            }
        }
    }
}