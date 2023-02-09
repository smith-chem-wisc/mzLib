using System;

namespace SimulatedData
{
    /// <summary>
    /// Extensions for the Random class. 
    /// </summary>
    internal static class RandomExtensions
    {
        /// <summary>
        /// Generates a random double. 
        /// </summary>
        /// <param name="random">Object of class Random.</param>
        /// <param name="low">The minimum value of the range.</param>
        /// <param name="high">The maximum value of the range.</param>
        /// <returns></returns>
        internal static double NextDouble(this Random random, double low, double high)
        {
            return random.NextDouble() * (high - low) + low;
        }
    }

    
}
