using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Oligo;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.Transcriptomics;

[TestFixture]
public static class BaseLossFragmentation
{
    [Test]
    public static void CanSuppressBaseLossIons_Enabled()
    {
        var unmodified = new OligoWithSetMods("GUACUG");
        var modified = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var modified2 = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]A[Biological:2'-O-Methyladenosine on A]UG");

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();
        var modifiedProducts2 = new List<Product>();

        FragmentationParams fragmentationParams = new RnaFragmentationParams
        {
             ModificationsCanSuppressBaseLossIons = true
        };

        unmodified.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts, fragmentationParams);
        modified.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts, fragmentationParams);
        modified2.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts2, fragmentationParams);

        Assert.That(unmodifiedProducts.Any(p => p.ProductType == ProductType.aBaseLoss));
        Assert.That(modifiedProducts.Any(p => p.ProductType == ProductType.aBaseLoss));
        Assert.That(modifiedProducts2.Any(p => p.ProductType == ProductType.aBaseLoss));

        Assert.That(unmodifiedProducts.Count, Is.EqualTo(modifiedProducts.Count + 1));
        Assert.That(unmodifiedProducts.Count, Is.EqualTo(modifiedProducts2.Count + 2));
    }

    [Test]
    public static void CanSuppressBaseLossIons_NotEnabled()
    {
        var unmodified = new OligoWithSetMods("GUACUG");
        var modified = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var modified2 = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]A[Biological:2'-O-Methyladenosine on A]UG");

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();
        var modifiedProducts2 = new List<Product>();

        FragmentationParams fragmentationParams = new RnaFragmentationParams
        {
            ModificationsCanSuppressBaseLossIons = false
        };

        unmodified.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts, fragmentationParams);
        modified.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts, fragmentationParams);
        modified2.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts2, fragmentationParams);

        Assert.That(unmodifiedProducts.Any(p => p.ProductType == ProductType.aBaseLoss));
        Assert.That(modifiedProducts.Any(p => p.ProductType == ProductType.aBaseLoss));
        Assert.That(modifiedProducts2.Any(p => p.ProductType == ProductType.aBaseLoss));

        Assert.That(unmodifiedProducts.Count, Is.EqualTo(modifiedProducts.Count));
        Assert.That(unmodifiedProducts.Count, Is.EqualTo(modifiedProducts2.Count));
    }

    [Test]
    public static void BaseLossDeltaMass_TwoOMethyl_IsCorrectForAllIons()
    {
        var unmodified = new OligoWithSetMods("GUACUG");
        var twoOMethyl = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;
        int modifiedResidue = 3;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        unmodified.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        twoOMethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        // Ensure all products except the base loss are identical
        var types = unmodifiedProducts.Select(p => p.ProductType).ToHashSet();
        foreach (var type in types)
        {
            if (type == ProductType.M)
            {
                continue;
            }
            var unmodifiedTypeProducts = unmodifiedProducts.Where(p => p.ProductType == type).ToList();
            var modifiedTypeProducts = modifiedProducts.Where(p => p.ProductType == type).ToList();

            Assert.That(unmodifiedTypeProducts.Count, Is.EqualTo(modifiedTypeProducts.Count));

            for (int i = 0; i < unmodifiedTypeProducts.Count; i++)
            {
                var unMod = unmodifiedTypeProducts[i];
                var mod = modifiedTypeProducts[i];
                var deltaMass = mod.NeutralMass - unMod.NeutralMass;

                Assert.That(unMod.FragmentNumber, Is.EqualTo(mod.FragmentNumber));
                Assert.That(unMod.ResiduePosition, Is.EqualTo(mod.ResiduePosition));
                Assert.That(unMod.Terminus, Is.EqualTo(mod.Terminus));

                if (type == ProductType.aBaseLoss)
                {
                    switch (i+1)
                    {
                        case 1: Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6)); break;
                        case 2: Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6)); break;
                        case 3: Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6)); break; // neutral loss will happen here
                        case 4: Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6)); break; 
                        case 5: Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6)); break;
                    }
                }
                else
                {
                    if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[type] == FragmentationTerminus.FivePrime)
                    {
                        if (mod.ResiduePosition < modifiedResidue)
                            Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6));
                        else
                            Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6));
                    }
                    else if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[type] == FragmentationTerminus.ThreePrime)
                    {
                        if (mod.ResiduePosition >= modifiedResidue)
                            Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6));
                        else
                            Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6));
                    }
                }
            }
        }
    }

    [Test]
    public static void BaseLossDeltaMass_TwoOMethyl_IsCorrectSimple()
    {
        var nSix = new OligoWithSetMods("GUACUG");
        var twoOMethyl = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        nSix.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        twoOMethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var twoOBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var nSixBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(twoOBaseLoss, Is.Not.Null);
        Assert.That(nSixBaseLoss, Is.Not.Null);

        var deltaMass = twoOBaseLoss.NeutralMass - nSixBaseLoss.NeutralMass;
        Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6));
    }

    [Test]
    public static void BaseLossDeltaMass_N6Methyl_IsCorrectSimple()
    {
        var unmodified = new OligoWithSetMods("GUACUG");
        var nSix = new OligoWithSetMods("GUA[N6-methyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        unmodified.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        nSix.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var twoOBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var nSixBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(twoOBaseLoss, Is.Not.Null);
        Assert.That(nSixBaseLoss, Is.Not.Null);

        var deltaMass = twoOBaseLoss.NeutralMass - nSixBaseLoss.NeutralMass;
        Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6));
    }

    [Test]
    public static void BaseLossDeltaMass_TwoOvsN6_IsCorrectSimple()
    {
        var nSix = new OligoWithSetMods("GUA[Biological:N6-methyladenosine on A]CUG");
        var twoOMethyl = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("CH2").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        nSix.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        twoOMethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var twoOBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var nSixBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(twoOBaseLoss, Is.Not.Null);
        Assert.That(nSixBaseLoss, Is.Not.Null);

        var deltaMass = twoOBaseLoss.NeutralMass - nSixBaseLoss.NeutralMass;
        Assert.That(deltaMass, Is.EqualTo(modMass).Within(1E-6));
    }

    [Test]
    public static void BaseLossDeltaMass_DiMethyl_IsCorrectSimple()
    {
        var unmod = new OligoWithSetMods("GUACUG");
        var dimethyl = new OligoWithSetMods("GUA[Biological:N6,2'-O-dimethyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("C2H4").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        unmod.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        dimethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var unmodBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var dimethylBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(unmodBaseLoss, Is.Not.Null);
        Assert.That(dimethylBaseLoss, Is.Not.Null);

        var deltaMass = Math.Abs(dimethylBaseLoss.NeutralMass - unmodBaseLoss.NeutralMass);
        Assert.That(deltaMass, Is.EqualTo(modMass / 2).Within(1E-6));
    }

    [Test]
    public static void BaseLossDeltaMass_TwoOvsDimethyl_IsCorrectSimple()
    {
        var twoOMethyl = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var dimethyl = new OligoWithSetMods("GUA[Biological:N6,2'-O-dimethyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("C2H4").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        twoOMethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        dimethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var twoOBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var dimethylBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(twoOBaseLoss, Is.Not.Null);
        Assert.That(dimethylBaseLoss, Is.Not.Null);

        var deltaMass = Math.Abs(dimethylBaseLoss.NeutralMass - twoOBaseLoss.NeutralMass);
        Assert.That(deltaMass, Is.EqualTo(0).Within(1E-6));
    }

    [Test]
    public static void BaseLossDeltaMass_N6vsDimethyl_IsCorrectSimple()
    {
        var nSix = new OligoWithSetMods("GUA[Biological:N6-methyladenosine on A]CUG");
        var dimethyl = new OligoWithSetMods("GUA[Biological:N6,2'-O-dimethyladenosine on A]CUG");
        var modMass = ChemicalFormula.ParseFormula("C2H4").MonoisotopicMass;

        var unmodifiedProducts = new List<Product>();
        var modifiedProducts = new List<Product>();

        nSix.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);
        dimethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, modifiedProducts);

        var nSixBaseLoss = modifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);
        var dimethylBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
            .FirstOrDefault(p => p.FragmentNumber == 3);

        Assert.That(nSixBaseLoss, Is.Not.Null);
        Assert.That(dimethylBaseLoss, Is.Not.Null);

        var deltaMass = Math.Abs(dimethylBaseLoss.NeutralMass - nSixBaseLoss.NeutralMass);
        Assert.That(deltaMass, Is.EqualTo(modMass / 2).Within(1E-6));
    }

    [Test]
    public static void dummy()
    {
        var unmodified = new OligoWithSetMods("GUACUG");
        var twoOMethyl = new OligoWithSetMods("GUA[Biological:2'-O-Methyladenosine on A]CUG");
        var nSix = new OligoWithSetMods("GUA[Biological:N6-methyladenosine on A]CUG");
        var dimethyl = new OligoWithSetMods("GUA[Biological:N6,2'-O-dimethyladenosine on A]CUG");

        var unmodifiedProducts = new List<Product>();
        unmodified.Fragment(DissociationType.CID, FragmentationTerminus.Both, unmodifiedProducts);

        var twoOMethylProducts = new List<Product>();
        twoOMethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, twoOMethylProducts);

        var nSixProducts = new List<Product>();
        nSix.Fragment(DissociationType.CID, FragmentationTerminus.Both, nSixProducts);

        var dimethylProducts = new List<Product>();
        dimethyl.Fragment(DissociationType.CID, FragmentationTerminus.Both, dimethylProducts);

        for (int i = 0; i < 5; i++)
        {
            var unmodBaseLoss = unmodifiedProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
                .FirstOrDefault(p => p.FragmentNumber == i + 1);
            var twoOBaseLoss = twoOMethylProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
                .FirstOrDefault(p => p.FragmentNumber == i + 1);
            var nSixBaseLoss = nSixProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
                .FirstOrDefault(p => p.FragmentNumber == i + 1);
            var dimethylBaseLoss = dimethylProducts.Where(p => p.ProductType == ProductType.aBaseLoss)
                .FirstOrDefault(p => p.FragmentNumber == i + 1);

            Console.WriteLine($"Fragment {i + 1}:");
            Console.WriteLine($"Unmodified: {unmodBaseLoss?.NeutralMass}");
            Console.WriteLine($"Two O-Methyl: {twoOBaseLoss?.NeutralMass}");
            Console.WriteLine($"N6 Methyl: {nSixBaseLoss?.NeutralMass}");
            Console.WriteLine($"Dimethyl: {dimethylBaseLoss?.NeutralMass}");
            Console.WriteLine();

        }


    }
}