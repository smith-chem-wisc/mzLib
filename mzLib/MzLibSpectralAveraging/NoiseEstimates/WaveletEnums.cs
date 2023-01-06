namespace MzLibSpectralAveraging;

/// <summary>
/// Options for the wavelets that are currently implemented. 
/// </summary>
/// <remarks>If expanding the supported wavelets, you will need to add the wavelet name here.</remarks>
public enum WaveletType
{
    Haar = 1, 
    Db4 = 2
}
/// <summary>
/// Options for boundary type in the modwt transfom. Only Reflection is supported right now. 
/// </summary>
public enum BoundaryType
{
    Reflection = 1
}