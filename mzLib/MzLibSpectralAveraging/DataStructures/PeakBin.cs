using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Xml.Serialization;

namespace MzLibSpectralAveraging;

/// <summary>
/// Class that stores and provides a means of performing operations such as rejection on a set of values in a pixel. 
/// </summary>
/// <remarks>Some development history: originally stored a double[] of the intensity values and tried to simultaneously keep track of
/// whether or not a value was rejected or not. This created some pretty complex logic and storage requirements. I realized that
/// what I was doing violated the single responsibility principle and refactored the class, breaking out the storage of m/z, intensity,
/// source spectra, and rejection status into the Pixel class. Now, this class only (1) performs operations on the pixel list and
/// (2) provides handles to the pixel list</remarks>
public class PeakBin
{
    public int Length => _peaks.Count;
    public int NonRejectedLength => _peaks.Count(i => i.Rejected == false); 
    public double MergedIntensityValue { get; private set; }

    /// <summary>
    /// Average of Unrejected Mz Values
    /// </summary>
    public double UnrejectedMzAverage => _peaks.Where(j => j.Rejected == false).Average(j => j.Mz);

    /// <summary>
    /// Intensity values of each pixel in stack
    /// </summary>
    public List<double> Intensities => _peaks.Select(i => i.Intensity).ToList();

    /// <summary>
    /// Mz values of each pixel in stack
    /// </summary>
    public List<double> Mzs => _peaks.Select(i => i.Mz).ToList();

    /// <summary>
    /// Intensity value of each unrejected pixel in stack
    /// </summary>
    public IEnumerable<double> UnrejectedIntensities => _peaks.Where(i => i.Rejected == false)
        .Select(i => i.Intensity);
   
    /// <summary>
    /// Mz Value of each unrejected pixel in stack
    /// </summary>
    public IEnumerable<double> UnrejectedMzs => _peaks.Where(i => i.Rejected == false)
        .Select(i => i.Mz);

    private List<Peak> _peaks { get; set; }

    public List<Peak> Peaks => _peaks;

    public PeakBin(IEnumerable<double> xArray, IEnumerable<double> yArray)
    {
        _peaks = xArray.Zip(yArray, (mz, its) => (mz,its))
            .Select((m,n) => new Peak(n, m.mz, m.its, rejected:false))
            .ToList();
        _peaks.Sort(new Peak.PixelComparer());
    }

    public PeakBin(IEnumerable<Peak> peaks)
    {
        _peaks = peaks.ToList();
        _peaks.Sort(new Peak.PixelComparer());
    }

    /// <summary>
    /// Sets the Rejection property of a Pixel object at the index specified to true. 
    /// </summary>
    /// <param name="index">The index of a pixel to be rejected.</param>
    public void Reject(int index)
    {
        _peaks[index].Rejected = true; 
    }

    /// <summary>
    /// Sets the Rejection property of all Pixel objects to true. 
    /// </summary>
    public void RejectAll()
    {
        _peaks.ForEach(p => p.Rejected = true);
    }

    /// <summary>
    /// Returns the Rejection value of the given pixel object
    /// </summary>
    /// <param name="index"></param>
    /// <returns>True is pixel is rejected. False if pixel is not rejected.</returns>
    public bool IsIndexRejected(int index)
    {
        return _peaks[index].Rejected;
    }

    /// <summary>
    /// Modifies the intensity property of a pixel value at a given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <param name="value"></param>
    internal void ModifyPixelIntensity(int index, double value)
    {
        _peaks[index].Intensity = value; 
    }

    /// <summary>
    /// Returns the intensity property of a pixel at the given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <returns></returns>
    public double GetIntensityAtIndex(int index)
    {
        return _peaks[index].Intensity; 
    }

    /// <summary>
    /// Returns the mz property of a pixel at the given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <returns></returns>
    public double GetMzAtIndex(int index)
    {
        return _peaks[index].Mz; 
    }

    /// <summary>
    /// Performs weighted average calculation givens a dictionary of weights.
    /// </summary>
    /// <param name="weightsDictionary">A dictionary of weights where the key is the source spectra and the value is the weight associated with that spectra.</param>
    /// <param name="originalTics"></param>
    /// 
    public void Average(IDictionary<int, double> weightsDictionary)
    {
        double numerator = 0;
        double denominator = 0;

        _peaks.Sort(new Peak.PixelComparer());

        foreach (var weight in weightsDictionary)
        {
            int index = _peaks.IndexOf(_peaks.First(i => i.SpectraId == weight.Key)); 
            if (_peaks[index].Rejected) continue;

            numerator += weight.Value * _peaks[index].Intensity;
            denominator += weight.Value; 
        }
        MergedIntensityValue = numerator / denominator;
    }

    public override string ToString()
    {
        return Mzs.Average() + " : " + Intensities.Average();
    }
}
