using System;
using System.Collections.Generic;
using Microsoft.Extensions.ObjectPool;

namespace MzLibUtil;


// Example Usage:
// var pool = new HashSetPool<int>();
// var hashSet = pool.Get();
// try {
// hashSet.Add(1);
// Do Work
// }
// finally {
// pool.Return(hashSet);
// }

/// <summary>
/// Provides a pool for <see cref="HashSet{T}"/> instances to reduce memory allocations.
/// This class uses the <see cref="ObjectPool{T}"/> from Microsoft.Extensions.ObjectPool
/// to manage the pooling of <see cref="HashSet{T}"/> objects.
/// </summary>
/// <typeparam name="T">The type of elements in the <see cref="HashSet{T}"/>.</typeparam>
/// <remarks>
/// This class is not thread-safe and should not be shared between threads.
/// This class should be pulled from outside a try finally loop and finally should return the HashSet to the pool to ensure proper pooling in the case of a caught exception
/// See example found in DigestionAgent.GetDigestionSiteIndices() for proper usage
/// </remarks>
public class HashSetPool<T>
{
    private readonly ObjectPool<HashSet<T>> _pool;

    /// <summary>
    /// Initializes a new instance of the <see cref="HashSetPool{T}"/> class.
    /// </summary>
    /// <param name="initialCapacity">Initial capacity for the pooled HashSet instances.</param>
    public HashSetPool(int initialCapacity = 16)
    {
        var policy = new HashSetPooledObjectPolicy<T>(initialCapacity);
        _pool = new DefaultObjectPool<HashSet<T>>(policy);
    }

    /// <summary>
    /// Retrieves a <see cref="HashSet{T}"/> instance from the pool.
    /// </summary>
    /// <returns>A <see cref="HashSet{T}"/> instance.</returns>
    public HashSet<T> Get() => _pool.Get();

    /// <summary>
    /// Returns a <see cref="HashSet{T}"/> instance back to the pool.
    /// </summary>
    /// <param name="hashSet">The <see cref="HashSet{T}"/> instance to return.</param>
    public void Return(HashSet<T> hashSet)
    {
        if (hashSet == null) throw new ArgumentNullException(nameof(hashSet));
        hashSet.Clear(); // Ensure the HashSet is clean before returning it to the pool
        _pool.Return(hashSet);
    }

    /// <summary>
    /// Defines the policy for creating and returning <see cref="HashSet{T}"/> instances to the pool.
    /// </summary>
    /// <typeparam name="TItem">The type of elements in the <see cref="HashSet{T}"/>.</typeparam>
    private class HashSetPooledObjectPolicy<TItem>(int initialCapacity) : PooledObjectPolicy<HashSet<TItem>>
    {
        /// <summary>
        /// Creates a new <see cref="HashSet{T}"/> instance with the specified initial capacity.
        /// </summary>
        /// <returns>A new <see cref="HashSet{T}"/> instance.</returns>
        public override HashSet<TItem> Create()
        {
            return new HashSet<TItem>(capacity: initialCapacity);
        }

        /// <summary>
        /// Returns a <see cref="HashSet{T}"/> instance to the pool after clearing it.
        /// </summary>
        /// <param name="obj">The <see cref="HashSet{T}"/> instance to return.</param>
        /// <returns>Always returns true.</returns>
        public override bool Return(HashSet<TItem> obj)
        {
            // Ensure the HashSet can be safely reused
            obj.Clear();
            return true;
        }
    }
}
