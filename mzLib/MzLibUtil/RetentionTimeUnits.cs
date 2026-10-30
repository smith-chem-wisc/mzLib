using System.Globalization;

namespace MzLibUtil
{
    /// <summary>
    /// An indexed retention time: the run-independent scale a spectral library stores (Biognosys/Prosit iRT). It is a
    /// different quantity from <see cref="RtMinutes"/>. There is deliberately no conversion between the two; a calibration
    /// fitted to one run is the only bridge.
    /// </summary>
    public readonly record struct Irt(double Value)
    {
        public override string ToString() => $"{Value.ToString(CultureInfo.InvariantCulture)} iRT";
    }

    /// <summary>A retention time in minutes within one run, as an MS scan reports it.</summary>
    public readonly record struct RtMinutes(double Value)
    {
        public override string ToString() => $"{Value.ToString(CultureInfo.InvariantCulture)} min";
    }
}
