﻿using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;

namespace SpectralAveraging;

public static class MassSpectrometryExtensions
{
    /// <summary>
    ///     Average an enumerable of MzSpectrum objects
    /// </summary>
    /// <param name="spectraToAverage">Spectra to average</param>
    /// <param name="parameters">Options for how to average spectra</param>
    /// <returns>Averaged MzSpectrum object</returns>
    public static MzSpectrum AverageSpectra(this IEnumerable<MzSpectrum> spectraToAverage,
        SpectralAveragingParameters parameters)
    {
        var toAverage = spectraToAverage as MzSpectrum[] ?? spectraToAverage.ToArray();
        var xArrays = toAverage.Select(p => p.XArray).ToArray();
        var yArrays = toAverage.Select(p => p.YArray.SubArray(0, p.YArray.Length)).ToArray();

        var xyJagged = SpectraAveraging.AverageSpectra(xArrays, yArrays, parameters);

        return new MzSpectrum(xyJagged[0], xyJagged[1], true);
    }

    /// <summary>
    /// Average an enumerable of MsDataScans objects
    /// </summary>
    /// <param name="scansToAverage">Scans to average</param>
    /// <param name="parameters">Options for how to average scans</param>
    /// <returns>Averaged MzSpectrum object</returns>
    public static MzSpectrum AverageSpectra(this IEnumerable<MsDataScan> scansToAverage,
        SpectralAveragingParameters parameters)
    {
        return scansToAverage.Select(p => p.MassSpectrum).AverageSpectra(parameters);
    }


    /// <summary>
    /// Normalize a group of MzSpectrum
    /// </summary>
    /// <param name="spectraToNormalize">spectra to normalize</param>
    /// <param name="type">normalization type</param>
    public static void NormalizeSpectra(this IEnumerable<MzSpectrum> spectraToNormalize, NormalizationType type)
    {
        var yArrays = spectraToNormalize.Select(p => p.YArray).ToArray();
        SpectraNormalization.NormalizeSpectra(yArrays, type);
    }

    /// <summary>
    /// Normalize a group of MsDataScans
    /// </summary>
    /// <param name="scansToNormalize">spectra to normalize</param>
    /// <param name="type">normalization type</param>
    public static void NormalizeSpectra(this IEnumerable<MsDataScan> scansToNormalize, NormalizationType type)
    {
        scansToNormalize.Select(p => p.MassSpectrum).NormalizeSpectra(type);
    }

    /// <summary>
    /// Normalize MzSpectrum
    /// </summary>
    /// <param name="spectrum">Spectrum to normalized</param>
    /// <param name="normalizationType">How To Normalize</param>
    public static void NormalizeSpectrum(this MzSpectrum spectrum, NormalizationType normalizationType)
    {
        var yArrays = new[] { spectrum.YArray };
        SpectraNormalization.NormalizeSpectra(yArrays, normalizationType);
    }

    /// <summary>
    /// Normalize a MsDataScan
    /// </summary>
    /// <param name="scan">Scan to normalize</param>
    /// <param name="normalizationType">How To Normalize</param>
    public static void NormalizeSpectrum(this MsDataScan scan, NormalizationType normalizationType)
    {
        scan.MassSpectrum.NormalizeSpectrum(normalizationType);
    }

    /// <summary>
    /// Converts an MzSpectrum y array to relative intensity
    /// </summary>
    /// <param name="spectrum">spectrum to convert</param>
    public static void ToRelativeIntensity(this MzSpectrum spectrum)
    {
        spectrum.NormalizeSpectrum(NormalizationType.RelativeIntensity);
    }

    /// <summary>
    /// Converts a MsDataScan to relative intensity
    /// </summary>
    /// <param name="scan">scan to convert</param>
    public static void ToRelativeIntensity(this MsDataScan scan)
    {
        scan.MassSpectrum.ToRelativeIntensity();
    }
}