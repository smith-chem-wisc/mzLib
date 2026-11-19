#nullable enable
namespace MassSpectrometry
{
    /// <summary>
    /// What a DIA precursor's quantity is in one run (QuantProject's DEF-DIA-KIND). Only <see cref="Value"/> carries a
    /// quantity; the other two are different kinds of "no value" and must not be read as 0.
    /// </summary>
    public enum PrecursorQuantityKind
    {
        /// <summary>
        /// A positive area from at least <see cref="FragmentQuantity.MinimumFragments"/> measured fragments within valid
        /// peak bounds, above the floor.
        /// </summary>
        Value,

        /// <summary>
        /// Valid peak bounds, but every measured fragment at or below the floor: the quantity is below what the run could
        /// show. It is censored, not zero and not missing.
        /// </summary>
        BelowFloor,

        /// <summary>
        /// Fewer than <see cref="FragmentQuantity.MinimumFragments"/> fragments could be measured, or the peak had no valid
        /// bounds. Nothing is known about the quantity.
        /// </summary>
        NotLooked,
    }

    /// <summary>A DIA precursor's quantity in one run: a value, a value below the floor, or nothing looked at.</summary>
    /// <param name="Value">The quantity when <paramref name="Kind"/> is <see cref="PrecursorQuantityKind.Value"/>, which is
    /// above <paramref name="Floor"/>; NaN otherwise, never 0.</param>
    /// <param name="Kind">What the quantity is.</param>
    /// <param name="Floor">The floor the fragments were judged against, in the quantity's units.</param>
    public readonly record struct PrecursorQuantity(double Value, PrecursorQuantityKind Kind, double Floor);
}
