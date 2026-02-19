using System;
using System.Collections.Generic;
using System.Threading;

namespace MzLibUtil.ObjectPools
{
    // Example Usage:
    // var list = SingletonDoublePool.Instance.Get();
    // try {
    //     list.Add(1.0);
    //     // Do Work
    // }
    // finally {
    //     SingletonDoublePool.Instance.Return(list);
    // }

    /// <summary>
    /// Provides a thread-safe singleton pool for <see cref="List{Double}"/> instances to reduce memory allocations.
    /// This class uses lazy initialization to ensure the pool is created only once on first access.
    /// </summary>
    /// <remarks>
    /// This class is thread-safe and can be shared between threads.
    /// Lists should be returned to the pool in a finally block to ensure proper pooling in the case of a caught exception.
    /// </remarks>
    public sealed class SingletonDoublePool
    {
        private static readonly Lazy<SingletonDoublePool> _instance = new Lazy<SingletonDoublePool>(() => new SingletonDoublePool());
        private readonly ListPool<double> _pool;

        /// <summary>
        /// Gets the singleton instance of the <see cref="SingletonDoublePool"/>.
        /// </summary>
        public static SingletonDoublePool Instance => _instance.Value;

        /// <summary>
        /// Initializes a new instance of the <see cref="SingletonDoublePool"/> class.
        /// Private constructor to enforce singleton pattern.
        /// </summary>
        private SingletonDoublePool()
        {
            ThreadPool.GetMaxThreads(out int maxThreads, out int completionPortThreads);
            _pool = new ListPool<double>(initialCapacity: maxThreads*2);
        }

        /// <summary>
        /// Retrieves a <see cref="List{Double}"/> instance from the pool.
        /// </summary>
        /// <returns>A <see cref="List{Double}"/> instance with initial capacity equal to twice the maximum number of threads in the thread pool</returns>
        public List<double> Get() => _pool.Get();

        /// <summary>
        /// Returns a <see cref="List{Double}"/> instance back to the pool.
        /// </summary>
        /// <param name="list">The <see cref="List{Double}"/> instance to return.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="list"/> is null.</exception>
        public void Return(List<double> list)
        {
            if (list == null) throw new ArgumentNullException(nameof(list));
            _pool.Return(list);
        }
    }
}
