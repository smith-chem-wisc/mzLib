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
public class PixelStack
{
    public int Length => _pixels.Count;
    public int NonRejectedLength => _pixels.Count(i => i.Rejected == false); 
    public double MergedIntensityValue { get; private set; }
    public double MergedMzValue => CalculateMzAverage(); 
    public List<double> Intensity => GetIntensities().ToList();
    public List<double> Mzs => GetMzValues().ToList(); 
    public IEnumerable<double> UnrejectedIntensities => GetUnrejectedIntValues();
    public IEnumerable<double> UnrejectedMzs => GetUnrejectedMzValues();  
    private List<Pixel> _pixels { get; set; }
    
    public PixelStack(IEnumerable<double> xArray, IEnumerable<double> yArray)
    {
        _pixels = xArray.Zip(yArray, (mz, its) => (mz,its))
            .Select((m,n) => new Pixel(n, m.mz, m.its, rejected:false))
            .ToList();
        _pixels.Sort(new Pixel.PixelComparer());
    }
    /// <summary>
    /// Get all intensities in a pixel stack, regardless of rejection status. 
    /// </summary>
    /// <returns></returns>
    private IEnumerable<double> GetIntensities()
    {
        return _pixels.Select(i => i.Intensity);
    }
    /// <summary>
    /// Gets all mz values in a pixel stack, regardless of rejection status. 
    /// </summary>
    /// <returns></returns>
    private IEnumerable<double> GetMzValues()
    {
        return _pixels.Select(i => i.Mz); 
    }
    /// <summary>
    /// Gets only intensities that are not rejected. 
    /// </summary>
    /// <returns>Unrejected intensity values.</returns>
    private IEnumerable<double> GetUnrejectedIntValues()
    {
        return _pixels.Where(i => i.Rejected == false)
            .Select(i => i.Intensity); 
    }
    /// <summary>
    /// Gets only mz values from pixels that aren't rejected. 
    /// </summary>
    /// <returns>Unrejected mz values.</returns>
    private IEnumerable<double> GetUnrejectedMzValues()
    {
        return _pixels.Where(i => i.Rejected == false)
            .Select(i => i.Mz);
    }
    /// <summary>
    /// Sets the Rejection property of a Pixel object at the index specified to true. 
    /// </summary>
    /// <param name="index">The index of a pixel to be rejected.</param>
    public void Reject(int index)
    {
        _pixels[index].Rejected = true; 
    }
    /// <summary>
    /// Returns the Rejection value of the given pixel object
    /// </summary>
    /// <param name="index"></param>
    /// <returns>True is pixel is rejected. False if pixel is not rejected.</returns>
    public bool IsIndexRejected(int index)
    {
        return _pixels[index].Rejected;
    }
    /// <summary>
    /// Modifies the intensity property of a pixel value at a given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <param name="value"></param>
    internal void ModifyPixelIntensity(int index, double value)
    {
        _pixels[index].Intensity = value; 
    }

    // not currently used. 
    //internal void ModifyPixelMz(int index, double value)
    //{
    //    _pixels[index].Intensity = value; 
    //}
    /// <summary>
    /// Returns the intensity property of a pixel at the given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <returns></returns>
    public double GetIntensityAtIndex(int index)
    {
        return _pixels[index].Intensity; 
    }

    /// <summary>
    /// Returns the mz property of a pixel at the given index. 
    /// </summary>
    /// <param name="index"></param>
    /// <returns></returns>
    public double GetMzAtIndex(int index)
    {
        return _pixels[index].Mz; 
    }
    /// <summary>
    /// Calculates the average mz value of only unrejected pixels. 
    /// </summary>
    /// <returns></returns>
    private double CalculateMzAverage()
    {
        return _pixels.Where(j => j.Rejected == false)
            .Average(j => j.Mz); 
    }
    /// <summary>
    /// Performs rejection on this pixel stack. 
    /// </summary>
    /// <param name="options">SpectralAveragingOptions object.</param>
    public void PerformRejection(SpectralAveragingOptions options)
    {
        OutlierRejection.RejectOutliers(this, options);
    }
    /// <summary>
    /// Performs weighted average calculation givens a dictionary of weights.
    /// </summary>
    /// <param name="weightsDictionary">A dictionary of weights where the key is the source spectra and the value is the weight associated with that spectra.</param>
    /// 
    public void Average(IDictionary<int, double> weightsDictionary)
    {
        double numerator = 0;
        double denominator = 0;

        _pixels.Sort(new Pixel.PixelComparer());

        foreach (var weight in weightsDictionary)
        {
            int index = _pixels.IndexOf(_pixels.Where(i => i.SpectraId == weight.Key).First()); 
            if (_pixels[index].Rejected == true) continue;

            numerator += weight.Value * _pixels[index].Intensity;
            denominator += weight.Value; 
        }
        MergedIntensityValue = numerator / denominator;
    }

    public override string ToString()
    {
        return Mzs.Average() + " : " + Intensity.Average();
    }
}
