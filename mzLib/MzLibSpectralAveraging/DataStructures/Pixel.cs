using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibSpectralAveraging; 
/// <summary>
/// Fundamental object that records an Mz and Intensity pair. Also maintains a record of the initial spectra that
/// the pair came from and whether or not this particular pixel has been rejected. 
/// </summary>
/// <remarks>This was the key object in maintaining the link between source spectra, data, and whether or not the data was rejected.
/// It was originally a record, but ultimately required me to use an object. However, I don't remember the reason why.</remarks>
internal class Pixel
{
    public int SpectraId;
    public double Intensity;
    public double Mz;
    public bool Rejected;
    public Pixel(int spectraId, double mz, double intensity, bool rejected)
    {
        SpectraId = spectraId;
        Intensity = intensity;
        Rejected = rejected;
        Mz = mz;
    }
    /// <summary>
    /// Custom comparer that facilitates searching and sorting lists of pixel objects. 
    /// </summary>
    internal class PixelComparer : IComparer<Pixel>
    {
        public int Compare(Pixel? x, Pixel? y)
        {
            return x.SpectraId.CompareTo(y.SpectraId);
        }
    }
}