using System;
using System.Collections.Generic;
using Microsoft.Extensions.ObjectPool;

namespace MzLibUtil;

// Example Usage:
// var pool = new ListPool<int>();
// var list = pool.Get();
// try {
// list.Add(1);
// Do Work
// }
// finally {
// pool.Return(list);
// }

/// <summary>
/// Provides a pool for <see cref="List{T}"/> instances to reduce memory allocations.
/// This class uses the <see cref="ObjectPool{T}"/> from Microsoft.Extensions.ObjectPool
/// to manage the pooling of <see cref="List{T}"/> objects.
/// </summary>
/// <typeparam name="T">The type of elements in the <see cref="List{T}"/>.</typeparam>
/// <remarks>
/// This class is not thread-safe and should not be shared between threads.
/// This class should be pulled from outside a try finally loop and finally should return the List to the pool to ensure proper pooling in the case of a caught exception.
/// </remarks>
public class ListPool<T>
{
    private readonly ObjectPool<List<T>> _pool;

    /// <summary>
    /// Initializes a new instance of the <see cref="ListPool{T}"/> class.
    /// </summary>
    /// <param name="initialCapacity">Initial capacity for the pooled HashSet instances.</param>
    public ListPool(int initialCapacity = 16)
    {
        var policy = new ListPooledObjectPolicy<T>(initialCapacity);
        var provider = new DefaultObjectPoolProvider { MaximumRetained = Environment.ProcessorCount * 2 };
        _pool = provider.Create(policy);
    }

    /// <summary>
    /// Retrieves a HashSet instance from the pool.
    /// </summary>
    /// <returns>A HashSet instance.</returns>
    public List<T> Get() => _pool.Get();

    /// <summary>
    /// Returns a HashSet instance back to the pool.
    /// </summary>
    /// <param name="list">The HashSet instance to return.</param>
    public void Return(List<T> list)
    {
        if (list == null) throw new ArgumentNullException(nameof(list));
        list.Clear(); // Ensure the HashSet is clean before returning it to the pool
        _pool.Return(list);
    }

    /// <summary>
    /// Policy for pooling List instances with a specified initial capacity.
    /// </summary>
    /// <typeparam name="TItem">The type of elements in the List.</typeparam>
    /// <param name="initialCapacity">The initial capacity for the pooled List instances.</param>
    private class ListPooledObjectPolicy<TItem>(int initialCapacity) : PooledObjectPolicy<List<TItem>>
    {
        private int InitialCapacity { get; } = initialCapacity;

        /// <summary>
        /// Creates a new List instance with the specified initial capacity.
        /// </summary>
        /// <returns>A new List instance.</returns>
        public override List<TItem> Create()
        {
            return new List<TItem>(capacity: InitialCapacity);
        }

        /// <summary>
        /// Resets the List instance to a clean state before returning it to the pool.
        /// </summary>
        /// <param name="obj">The List instance to reset and return.</param>
        /// <returns>True if the List instance can be returned to the pool; otherwise, false.</returns>
        public override bool Return(List<TItem> obj)
        {
            // Ensure the List can be safely reused
            obj.Clear();
            return true;
        }
    }
}