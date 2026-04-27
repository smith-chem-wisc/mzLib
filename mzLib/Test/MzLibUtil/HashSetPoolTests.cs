using MzLibUtil;
using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;

namespace Test.MzLibUtil;

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
}
