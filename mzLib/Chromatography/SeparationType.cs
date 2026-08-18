namespace Chromatography;

/// <summary>
/// Types of chromatographic or electrophoretic separations supported for retention time prediction
/// </summary>
public enum SeparationType
{
    /// <summary>
    /// Reverse phase HPLC (High Performance Liquid Chromatography)
    /// Used by: SSRCalc3, Chronologer
    /// </summary>
    HPLC,
    
    /// <summary>
    /// CZE (Capillary Zone Electrophoresis)
    /// Used by: CZE predictor
    /// </summary>
    CZE
}
