using MzLibUtil;
using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;

namespace Test.Util;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ListPoolTests
{
    [Test]
    public void ListPool_Get_ReturnsListWithInitialCapacity()
    {
        int initialCapacity = 16;
        var listPool = new ListPool<int>(initialCapacity);

        var list = listPool.Get();

        Assert.That(list, Is.Not.Null);
        Assert.That(list.Capacity, Is.EqualTo(initialCapacity));
    }

    [Test]
    public void ListPool_Return_ClearsListBeforeReturningToPool()
    {
        var listPool = new ListPool<int>();
        var list = listPool.Get();
        list.Add(1);
        list.Add(2);

        listPool.Return(list);
        var returnedList = listPool.Get();

        Assert.That(returnedList, Is.Not.Null);
        Assert.That(returnedList, Is.Empty);
    }

    [Test]
    public void ListPool_Return_ThrowsArgumentNullException_WhenListIsNull()
    {
        var listPool = new ListPool<int>();

        Assert.That(() => listPool.Return(null), Throws.ArgumentNullException);
    }

    /// <summary>
    /// Nothing asserted that the pool actually pools. The policy's Return can refuse an instance, in
    /// which case Get quietly allocates a fresh one -- observationally identical unless the test looks
    /// at instance identity.
    /// </summary>
    [Test]
    public void ListPool_Return_MakesTheSameInstanceAvailableAgain()
    {
        var listPool = new ListPool<int>();
        var list = listPool.Get();
        list.Add(1);

        listPool.Return(list);

        Assert.That(listPool.Get(), Is.SameAs(list));
    }

}
