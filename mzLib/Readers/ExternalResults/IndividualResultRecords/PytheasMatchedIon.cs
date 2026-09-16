using System;
using System.Globalization;
using System.Text.RegularExpressions;
using MzLibUtil;

namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// One normalized MS2 match printed in a Pytheas match record, e.g.
/// "892.167114(0.8ppm)[680]:892.166417[M-P(-1)]". The token reports a measured peak m/z and its
/// matched theoretical fragment: the measured/offset/intensity triple before the ':' and the
/// theoretical m/z with its CID-ion annotation after it. Parsed here so callers get typed numbers
/// instead of the raw string; the original token is kept so the record can round-trip on write.
/// </summary>
public class PytheasMatchedIon
{
    private static readonly Regex TokenRegex = new(
        @"^([\d.]+)\(\s*(-?[\d.]+)\s*ppm\)\[([\d.]+)]:([\d.]+)\[(.+)]$",
        RegexOptions.Compiled);

    private static readonly Regex ChargeRegex = new(@"\(([+-]\d+)\)$", RegexOptions.Compiled);

    /// <summary>The measured m/z of the fragment peak (left of the ':').</summary>
    public double MeasuredMz { get; }

    /// <summary>offset: measured minus theoretical m/z, in ppm.</summary>
    public double OffsetPpm { get; }

    /// <summary>norm_intensity: normalized intensity of the fragment peak.</summary>
    public double NormalizedIntensity { get; }

    /// <summary>The theoretical m/z of the matched fragment (right of the ':').</summary>
    public double TheoreticalMz { get; }

    /// <summary>The CID-ion annotation, e.g. "M-P(-1)", "y1(-1)", "c4(-1)", "M-H2O-P(-1)".</summary>
    public string IonAnnotation { get; }

    /// <summary>The charge parsed from the annotation's trailing "(+/-n)" group.</summary>
    public int Charge { get; }

    /// <summary>The original token, kept verbatim for round-trip writing.</summary>
    public string RawMatch { get; }

    public PytheasMatchedIon(double measuredMz, double offsetPpm, double normalizedIntensity,
        double theoreticalMz, string ionAnnotation, int charge, string rawMatch)
    {
        MeasuredMz = measuredMz;
        OffsetPpm = offsetPpm;
        NormalizedIntensity = normalizedIntensity;
        TheoreticalMz = theoreticalMz;
        IonAnnotation = ionAnnotation;
        Charge = charge;
        RawMatch = rawMatch;
    }

    public static PytheasMatchedIon Parse(string token)
    {
        var match = TokenRegex.Match(token);
        if (!match.Success)
            throw new MzLibException($"Could not parse Pytheas MS2 match token: {token}");

        string annotation = match.Groups[5].Value;
        var chargeMatch = ChargeRegex.Match(annotation);
        int charge = chargeMatch.Success
            ? int.Parse(chargeMatch.Groups[1].Value, CultureInfo.InvariantCulture)
            : 0;

        return new PytheasMatchedIon(
            measuredMz: ParseDouble(match.Groups[1].Value),
            offsetPpm: ParseDouble(match.Groups[2].Value),
            normalizedIntensity: ParseDouble(match.Groups[3].Value),
            theoreticalMz: ParseDouble(match.Groups[4].Value),
            ionAnnotation: annotation,
            charge: charge,
            rawMatch: token);
    }

    private static double ParseDouble(string value) =>
        double.Parse(value, NumberStyles.Float, CultureInfo.InvariantCulture);
}