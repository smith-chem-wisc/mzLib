using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    /// <summary>
    /// Contract for an external deconvolution result used to pair with raw MS2 scans
    /// and (optionally) construct external MS2 spectra for database searching.
    /// Implementations are typically produced by mzLib Readers that parse
    /// third-party deconvolution outputs.
    /// </summary>
    /// <remarks>
    /// Required fields:
    /// <list type="bullet">
    /// <item><description><see cref="Mz"/>: precursor m/z (double)</description></item>
    /// <item><description><see cref="Charge"/>: precursor charge state (int)</description></item>
    /// <item><description><see cref="RetentionTimeStart"/>: RT window start, minutes (double)</description></item>
    /// <item><description><see cref="RetentionTimeEnd"/>: RT window end, minutes (double)</description></item>
    /// <item><description><see cref="NumberOfIsotopes"/>: Number of isotopes observed, nullable int</description></item>
    /// </list>
    /// Pairing is performed by checking if an MS2 scan’s RT apex lies within
    /// [<see cref="RetentionTimeStart"/>, <see cref="RetentionTimeEnd"/>].
    /// Implementations should ensure values are valid (no NaN; charge &gt; 0; start ≤ end).
    /// </remarks>
    public interface ISingleChargeMs1Feature
    {
        double Mz { get; }
        int Charge { get; }
        double RetentionTimeStart { get; }
        double RetentionTimeEnd { get; }
        double Intensity { get; }
        int? NumberOfIsotopes { get; }
    }
}
