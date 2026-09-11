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
    /// Every model's session deadline must cover the work it is bounding, at every request size.
    /// </summary>
    /// <remarks>
    /// This is the offline guard that lets <see cref="KoinaTests.KoinaLiveTestFixture"/> skip a test
    /// whose Koina call was cut off. A call that died because Koina never answered and a call that
    /// died because our own deadline was too short abort at the same line and read identically, so
    /// that fixture cannot tell them apart and skips both. Which is only safe while the second cannot
    /// happen -- and that is what this asserts.
    ///
    /// Arithmetic rather than network on purpose. A batching regression has to fail a deterministic
    /// test in the required job, not depend on a live test in the non-blocking one to notice it; a
    /// live test that goes quietly Skipped is precisely the outcome that got "out of memory" removed
    /// from KoinaServiceException.ServiceFaultMarkers.
    ///
    /// The bound asserted is the estimate's own definition -- twice the benchmarked time for every
    /// batch, plus the throttling delay between chunks -- so raising a batch size, lowering a
    /// benchmark, or rounding the deadline down all fail here rather than in CI a week later.
    /// </remarks>
    [Test]
    [TestCaseSource(nameof(ConcreteKoinaModelTypes))]
    public static void EveryConcreteModel_SessionDeadlineCoversItsOwnBenchmarkedWork(Type modelType)
    {
        var model = InstantiateWithDefaults(modelType);

        int throttle = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.ThrottlingDelayInMilliseconds))!.GetValue(model)!;
        int benchmark = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.BenchmarkedTimeForOneMaxBatchSizeInMilliseconds))!.GetValue(model)!;
        int maxBatchesPerRequest = (int)modelType.GetProperty(nameof(KoinaModelBase<int, int>.MaxNumberOfBatchesPerRequest))!.GetValue(model)!;
        var sessionDeadline = modelType.GetMethod(nameof(KoinaModelBase<int, int>.SessionDeadline))!;

        // One batch is the shape every live Koina test takes (a handful of peptides), and the shape
        // the 2026-09 AlphaPeptDeep_ccs_generic stall was observed on. The large counts are where a rounding or tuning
        // regression would actually bite.
        foreach (int batchCount in new[] { 1, 2, 17, 250, 1_000, 25_000 })
        {
            int chunkCount = (int)Math.Ceiling(batchCount / (double)maxBatchesPerRequest);
            var deadline = (TimeSpan)sessionDeadline.Invoke(model, new object[] { batchCount, chunkCount })!;

            double estimatedMilliseconds = batchCount * 2.0 * benchmark + (double)throttle * chunkCount;
            Assert.That(deadline.TotalMilliseconds, Is.GreaterThanOrEqualTo(estimatedMilliseconds),
                $"{modelType.Name}: the deadline for {batchCount} batch(es) in {chunkCount} chunk(s) is "
                + "shorter than the work it bounds, so a healthy run can be cut off -- and the live "
                + "fixture would report that as a Koina outage.");

            Assert.That(deadline, Is.GreaterThanOrEqualTo(TimeSpan.FromMinutes(1)),
                $"{modelType.Name}: the one-minute floor is the only thing standing between a small "
                + "request and a deadline measured in milliseconds.");
        }
    }

    /// <summary>
    /// Live (NETWORK) check: every model's ModelName must name a real endpoint on the Koina server.
    /// A model sends its ModelName verbatim as the request URL path, so a single typo (for example an
    /// extra underscore) silently turns every prediction for that model into a 404 at runtime. This
    /// case catches that against the real registry instead of trusting the spelling in the code.
    ///
    /// Tagged [Category("ExternalService")] so it runs in the dedicated non-blocking job rather than the
    /// required one, and routed through <see cref="ExternalServiceTestHelper.RunAsync"/> so an unreachable
    /// Koina server reports Skipped with a reason instead of the false failure this test used to produce.
    /// A wrong model name still fails, because the assertion below is not an availability problem.
    /// [TestCaseSource] runs it once per discovered model, so a wrong name fails as its own case.
    ///
    /// The rest of this fixture is offline and deliberately stays in the required job.
    /// </summary>
    [Test]
    [Category("ExternalService")]
    [Category("Koina")]
    [TestCaseSource(nameof(ConcreteKoinaModelTypes))]
    public static Task EveryConcreteModel_ModelNameResolvesToRegisteredKoinaEndpoint(Type modelType) =>
        ExternalServiceTestHelper.RunAsync("Koina", async () =>
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

        // Separate "Koina is unwell" from "this name is wrong" before asserting: 408/429/5xx are the
        // server's problem and skip, where the 4xx it returns for an unknown name must still fail.
        ExternalServiceTestHelper.ThrowIfUnavailable(response);

        // A success status means the name is registered. On failure the message names the offending
        // class, its ModelName, and the HTTP status so the fix is obvious — compare against the docs
        // and keep ModelName, the doc-URL comment, and the PR table in agreement.
        Assert.That(response.IsSuccessStatusCode, Is.True,
            $"{modelType.Name}.ModelName '{modelName}' did not resolve at Koina " +
            $"({(int)response.StatusCode} {response.ReasonPhrase}). Verify the exact identifier against " +
            "https://koina.wilhelmlab.org/docs.");
    });
}
