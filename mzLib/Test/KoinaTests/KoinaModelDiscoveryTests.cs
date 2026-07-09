using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http;
using System.Reflection;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;

namespace Test.KoinaTests;

/// <summary>
/// Reflects over the PredictionClients assembly to discover every concrete Koina model
/// and asserts that each one is constructible with its default ctor and exposes
/// well-formed required metadata. Catches "forgot to declare ModelName" or
/// "MaxBatchSize defaulted to 0" mistakes for any future concrete model without
/// requiring per-model boilerplate.
///
/// The offline cases validate that metadata is well-formed. One additional case
/// (<see cref="EveryConcreteModel_ModelNameResolvesToRegisteredKoinaEndpoint"/>) goes a step
/// further and verifies, over the network, that each ModelName actually exists on the Koina
/// server — it is tagged [Category("Koina")] so it only runs when network tests are enabled.
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

    /// <summary>
    /// Live (NETWORK) check: every model's ModelName must name a real endpoint on the Koina server.
    /// A model sends its ModelName verbatim as the request URL path, so a single typo (for example an
    /// extra underscore) silently turns every prediction for that model into a 404 at runtime. This
    /// case catches that against the real registry instead of trusting the spelling in the code.
    ///
    /// Tagged [Category("Koina")] so it only runs when network tests are enabled; it will report a
    /// false failure if the Koina server itself is unreachable. [TestCaseSource] runs it once per
    /// discovered model, so a wrong name fails as its own clearly-labelled case.
    /// </summary>
    [Test]
    [Category("Koina")]
    [TestCaseSource(nameof(ConcreteKoinaModelTypes))]
    public static async Task EveryConcreteModel_ModelNameResolvesToRegisteredKoinaEndpoint(Type modelType)
    {
        // ConcreteKoinaModelTypes() (above) discovered this Type by reflection — we only have a Type,
        // not a statically-typed reference. Build an instance with its default ctor, then read the
        // ModelName property off that instance reflectively. The "!" marks the values as known non-null.
        var model = InstantiateWithDefaults(modelType);
        string modelName = (string)modelType.GetProperty(nameof(KoinaModelBase<int, int>.ModelName))!.GetValue(model)!;

        // A short-lived client just for this probe — NOT the shared production HTTP.Client. "using"
        // disposes it when the method returns; 30s is plenty for a tiny metadata GET.
        using var client = new HttpClient { Timeout = TimeSpan.FromSeconds(30) };

        // HTTP.ModelsURL is "https://koina.wilhelmlab.org:443/v2/models/", so this requests
        // ".../v2/models/{modelName}", the model's metadata entry. The server answers 200 for a
        // registered model and a 4xx (e.g. 400) for a name it does not recognize.
        using var response = await client.GetAsync($"{HTTP.ModelsURL}{modelName}");

        // A success status means the name is registered. On failure the message names the offending
        // class, its ModelName, and the HTTP status so the fix is obvious — compare against the docs
        // and keep ModelName, the doc-URL comment, and the PR table in agreement.
        Assert.That(response.IsSuccessStatusCode, Is.True,
            $"{modelType.Name}.ModelName '{modelName}' did not resolve at Koina " +
            $"({(int)response.StatusCode} {response.ReasonPhrase}). Verify the exact identifier against " +
            "https://koina.wilhelmlab.org/docs.");
    }
}
