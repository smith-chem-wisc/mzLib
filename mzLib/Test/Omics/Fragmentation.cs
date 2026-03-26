using System;
using System.Linq;
using NUnit.Framework;
using Omics.Fragmentation;

namespace Test.Omics;
[TestFixture]
public static class Fragmentation
{
    [Test]
    public static void Extensions_IsBaseLoss()
    {
        foreach (var productType in Enum.GetValues<ProductType>())
        {
            if (productType.ToString().Contains("BaseLoss"))
                Assert.That(productType.IsBaseLoss(), Is.True, $"{productType} should be identified as a base loss ion.");
            else
                Assert.That(productType.IsBaseLoss(), Is.False, $"{productType} should not be identified as a base loss ion.");
        }   
    }

    [Test]
    public static void Extensions_IsWaterLoss()
    {
        foreach (var productType in Enum.GetValues<ProductType>())
        {
            if (productType.ToString().Contains("WaterLoss"))
                Assert.That(productType.IsWaterLoss(), Is.True, $"{productType} should be identified as a water loss ion.");
            else
                Assert.That(productType.IsWaterLoss(), Is.False, $"{productType} should not be identified as a water loss ion.");
        }
    }

    [Test]
    public static void Extensions_GetAssociatedTypes()
    {
        // Test that a, a*, a°, a-H2O, and a-BaseLoss are all associated with each other
        foreach (var type in Enum.GetValues<ProductType>())
        {
            var letterStarting = Enum.GetValues<ProductType>().Where(p => p.ToString().StartsWith(type.ToString()[0]) && p != type).ToArray();

            var associatedTypes = type.GetFragmentFamilyMembers();

            Assert.That(associatedTypes, Is.EquivalentTo(letterStarting), $"Associated types for {type} should be {string.Join(", ", letterStarting)}, but got {string.Join(", ", associatedTypes)}");
        }
    }
}
