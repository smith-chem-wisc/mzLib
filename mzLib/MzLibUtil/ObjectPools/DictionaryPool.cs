using System;
using System.Collections.Generic;
using Microsoft.Extensions.ObjectPool;

namespace MzLibUtil;

// Used to pool HashSet instances to reduce memory allocations
public class DictionaryPool<TKey, TValue> where TKey : notnull
{
    private readonly ObjectPool<Dictionary<TKey, TValue>> _pool;

    /// <summary>
    /// Initializes a new instance of the <see cref="DictionaryPool{TKey, TValue}"/> class.
    /// </summary>
    /// <param name="initialCapacity">Initial capacity for the pooled Dictionary instances.</param>
    public DictionaryPool(int initialCapacity = 16)
    {
        var policy = new DictionaryPooledObjectPolicy<TKey, TValue>(initialCapacity);
        var provider = new DefaultObjectPoolProvider { MaximumRetained = Environment.ProcessorCount * 2 };
        _pool = provider.Create(policy);
    }

    /// <summary>
    /// Retrieves a Dictionary instance from the pool.
    /// </summary>
    /// <returns>A Dictionary instance.</returns>
    public Dictionary<TKey, TValue> Get() => _pool.Get();

    /// <summary>
    /// Returns a Dictionary instance back to the pool.
    /// </summary>
    /// <param name="dictionary">The Dictionary instance to return.</param>
    public void Return(Dictionary<TKey, TValue> dictionary)
    {
        if (dictionary == null) throw new ArgumentNullException(nameof(dictionary));
        dictionary.Clear(); // Ensure the Dictionary is clean before returning it to the pool
        _pool.Return(dictionary);
    }

    private class DictionaryPooledObjectPolicy<TKeyItem, TValueItem>(int initialCapacity)
        : PooledObjectPolicy<Dictionary<TKeyItem, TValueItem>>
        where TKeyItem : notnull
    {
        private int InitialCapacity { get; } = initialCapacity;

        public override Dictionary<TKeyItem, TValueItem> Create()
        {
            return new Dictionary<TKeyItem, TValueItem>(capacity: InitialCapacity);
        }

        public override bool Return(Dictionary<TKeyItem, TValueItem> obj)
        {
            // Ensure the Dictionary can be safely reused
            obj.Clear();
            return true;
        }
    }
}