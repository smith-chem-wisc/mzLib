using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;

namespace Test.Util;

[TestFixture]
[ExcludeFromCodeCoverage]
public class DictionaryPoolTests
{
    [Test]
    public void Get_ReturnsDictionaryInstance()
    {
        var dictionaryPool = new DictionaryPool<string, int>();
        var dictionary = dictionaryPool.Get();
        Assert.That(dictionary, Is.Not.Null);
        Assert.That(dictionary, Is.InstanceOf<Dictionary<string, int>>());
    }

    [Test]
    public void Return_ClearsAndReturnsDictionaryToPool()
    {
        var dictionaryPool = new DictionaryPool<string, int>();
        var dictionary = dictionaryPool.Get();
        dictionary["key"] = 42;

        dictionaryPool.Return(dictionary);

        Assert.That(dictionary.Count, Is.EqualTo(0));
    }

    [Test]
    public void Return_ThrowsArgumentNullException_WhenDictionaryIsNull()
    {
        var dictionaryPool = new DictionaryPool<string, int>();
        Assert.Throws<ArgumentNullException>(() => dictionaryPool.Return(null));
    }
}
