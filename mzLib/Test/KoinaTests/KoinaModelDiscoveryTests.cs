using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;

namespace Test.KoinaTests;

/// <summary>
/// Reflects over the PredictionClients assembly to discover every concrete Koina model
/// and asserts that each one is constructible with its default ctor and exposes
/// well-formed required metadata. Catches "forgot to declare ModelName" or
/// "MaxBatchSize defaulted to 0" mistakes for any future concrete model without
/// requiring per-model boilerplate.
/// </summary>
[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class KoinaModelDiscoveryTests
{
    private static IEnumerable<Type> ConcreteKoinaModelTypes()
    {
        var assembly = typeof(FragmentIntensityModel).Assembly;
        foreach (var type in assembly.GetTypes())
        {
            if (type.IsAbstract || type.IsInterface)
            {
                continue;
            }

            // Walk up the inheritance chain to find KoinaModelBase<,>.
            for (var t = type.BaseType; t != null; t = t.BaseType)
            {
                if (t.IsGenericType && t.GetGenericTypeDefinition() == typeof(KoinaModelBase<,>))
                {
                    yield return type;
                    break;
                }
            }
        }
    }

    private static object InstantiateWithDefaults(Type modelType)
    {
        var ctor = modelType.GetConstructors()
            .FirstOrDefault(c => c.GetParameters().All(p => p.HasDefaultValue));
        if (ctor == null)
        {
            throw new InvalidOperationException(
                $"{modelType.FullName} has no constructor where all parameters have default values.");
        }
        var args = ctor.GetParameters().Select(p => p.DefaultValue).ToArray();
        return ctor.Invoke(args);
    }

    [Test]
    public static void EveryConcreteModel_HasNoRequiredConstructorArguments()
    {
        var missing = ConcreteKoinaModelTypes()
            .Where(t => t.GetConstructors().All(c => c.GetParameters().Any(p => !p.HasDefaultValue)))
            .Select(t => t.FullName)
            .ToList();

        Assert.That(missing, Is.Empty,
            "Every concrete Koina model must be constructible with no required arguments (use default parameter values).");
    }

    [Test]
    [TestCaseSource(nameof(ConcreteKoinaModelTypes))]
    public static void EveryConcreteModel_HasWellFormedMetadata(Type modelType)
    {
        object model;
        try
        {
            model = InstantiateWithDefaults(modelType);
        }
        catch (Exception ex)
        {
            Assert.Fail($"Cannot instantiate {modelType.FullName}: {ex.GetBaseException().Message}");
            return;
        }

        string modelName = (string)modelType.GetProperty(nameof(KoinaModelBase<int, int>.ModelName))!.GetValue(model)!;
        Assert.That(modelName, Is.Not.Null.And.Not.Empty, $"{modelType.Name}.ModelName must be a non-empty string.");

        int maxBatch = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.MaxBatchSize))!.GetValue(model)!;
        Assert.That(maxBatch, Is.GreaterThan(0), $"{modelType.Name}.MaxBatchSize must be positive.");

        int maxBatchesPerReq = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.MaxNumberOfBatchesPerRequest))!.GetValue(model)!;
        Assert.That(maxBatchesPerReq, Is.GreaterThan(0), $"{modelType.Name}.MaxNumberOfBatchesPerRequest must be positive (init from default ctor).");

        int throttle = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.ThrottlingDelayInMilliseconds))!.GetValue(model)!;
        Assert.That(throttle, Is.GreaterThanOrEqualTo(0), $"{modelType.Name}.ThrottlingDelayInMilliseconds must be non-negative.");

        int benchmark = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.BenchmarkedTimeForOneMaxBatchSizeInMilliseconds))!.GetValue(model)!;
        Assert.That(benchmark, Is.GreaterThan(0), $"{modelType.Name}.BenchmarkedTimeForOneMaxBatchSizeInMilliseconds must be positive.");

        int minLen = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.MinPeptideLength))!.GetValue(model)!;
        int maxLen = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.MaxPeptideLength))!.GetValue(model)!;
        Assert.That(minLen, Is.GreaterThanOrEqualTo(1), $"{modelType.Name}.MinPeptideLength must be >= 1.");
        Assert.That(maxLen, Is.GreaterThanOrEqualTo(minLen), $"{modelType.Name}.MaxPeptideLength must be >= MinPeptideLength.");
    }
}
