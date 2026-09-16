using MzLibUtil;
using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;

namespace Test.Util;

[TestFixture]
[ExcludeFromCodeCoverage]
public class HashSetPoolTests
{
    [Test]
    public void Get_ReturnsHashSetInstance()
    {
        var pool = new HashSetPool<int>();
        var hashSet = pool.Get();
        Assert.That(hashSet, Is.Not.Null);
        pool.Return(hashSet);
    }

    [Test]
    public void Return_ClearsHashSetBeforeReturningToPool()
    {
        var pool = new HashSetPool<int>();
        var hashSet = pool.Get();
        hashSet.Add(1);
        pool.Return(hashSet);
        Assert.That(hashSet.Count, Is.EqualTo(0));
    }

    [Test]
    public void Return_ThrowsArgumentNullException_WhenHashSetIsNull()
    {
        var pool = new HashSetPool<int>();
        Assert.Throws<ArgumentNullException>(() => pool.Return(null));
    }

    /// <summary>
    /// Nothing asserted that the pool actually pools. The policy's Return can refuse an instance, in
    /// which case Get quietly allocates a fresh one -- observationally identical unless the test looks
    /// at instance identity.
    /// </summary>
    [Test]
    public void HashSetPool_Return_MakesTheSameInstanceAvailableAgain()
    {
        var hashSetPool = new HashSetPool<int>();
        var hashSet = hashSetPool.Get();
        hashSet.Add(1);

        hashSetPool.Return(hashSet);

        Assert.That(hashSetPool.Get(), Is.SameAs(hashSet));
    }

}
