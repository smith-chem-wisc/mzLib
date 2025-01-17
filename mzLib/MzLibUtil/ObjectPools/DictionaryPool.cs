using System;
using System.Collections.Generic;
using Microsoft.Extensions.ObjectPool;

namespace MzLibUtil;

// Example Usage:
// var pool = new DictionaryPool<int, int>();
// var dictionary = pool.Get();
// try {
// dictionary.Add(1,1);
// Do Work
// }
// finally {
// pool.Return(dictionary);
// }

/// <summary>
/// Provides a pool for <see cref="Dictionary{TKey, TValue}"/> instances to reduce memory allocations.
/// This class uses the <see cref="ObjectPool{T}"/> from Microsoft.Extensions.ObjectPool
/// to manage the pooling of <see cref="Dictionary{TKey, TValue}"/> objects.
/// </summary>
/// <typeparam name="TKey">The type of keys in the <see cref="Dictionary{TKey, TValue}"/>.</typeparam>
/// <typeparam name="TValue">The type of values in the <see cref="Dictionary{TKey, TValue}"/>.</typeparam>
/// <remarks>
/// This class is not thread-safe and should not be shared between threads.
/// This class should be pulled from outside a try finally loop and finally should return the Dictionary to the pool to ensure proper pooling in the case of a caught exception.
/// </remarks>
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

    /// <summary>
    /// Policy for pooling Dictionary instances with a specified initial capacity.
    /// </summary>
    /// <typeparam name="TKeyItem">The type of keys in the Dictionary.</typeparam>
    /// <typeparam name="TValueItem">The type of values in the Dictionary.</typeparam>
    /// <param name="initialCapacity">The initial capacity for the pooled Dictionary instances.</param>
    private class DictionaryPooledObjectPolicy<TKeyItem, TValueItem>(int initialCapacity)
        : PooledObjectPolicy<Dictionary<TKeyItem, TValueItem>>
        where TKeyItem : notnull
    {
        private int InitialCapacity { get; } = initialCapacity;

        /// <summary>
        /// Creates a new Dictionary instance with the specified initial capacity.
        /// </summary>
        /// <returns>A new Dictionary instance.</returns>
        public override Dictionary<TKeyItem, TValueItem> Create()
        {
            return new Dictionary<TKeyItem, TValueItem>(capacity: InitialCapacity);
        }

        /// <summary>
        /// Returns a Dictionary instance to the pool after clearing it.
        /// </summary>
        /// <param name="obj">The Dictionary instance to return.</param>
        /// <returns>True if the Dictionary instance can be reused; otherwise, false.</returns>
        public override bool Return(Dictionary<TKeyItem, TValueItem> obj)
        {
            // Ensure the Dictionary can be safely reused
            obj.Clear();
            return true;
        }
    }
}
