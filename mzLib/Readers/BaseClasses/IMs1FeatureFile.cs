using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    /// <summary>
    /// Contract for a reader or adapter that can produce a sequence of
    /// <see cref="ISingleChargeMs1Feature"/> instances from an external
    /// deconvolution results file.
    /// </summary>
    /// <remarks>
    /// Implementations are responsible for:
    /// <list type="bullet">
    /// <item>
    ///     <description>
    ///     Parsing the file format (e.g., CSV, vendor export, mzML with decon results)
    ///     into one or more <see cref="ISingleChargeMs1Feature"/> objects.
    ///     </description>
    /// </item>
    /// <item>
    ///     <description>
    ///     Mapping source file fields to the required interface properties:
    ///     <see cref="ISingleChargeMs1Feature.MZ"/>,
    ///     <see cref="ISingleChargeMs1Feature.Charge"/>,
    ///     <see cref="ISingleChargeMs1Feature.RetentionTimeStart"/>,
    ///     <see cref="ISingleChargeMs1Feature.RetentionTimeEnd"/>,
    ///     <see cref="ISingleChargeMs1Feature.NumberOfIsotopes"/>.
    ///     </description>
    /// </item>
    /// </list>
    /// This interface provides a consistent entry point for the pairing logic
    /// to consume external results from various deconvolution tools without
    /// caring about the underlying file format.
    /// </remarks>
    public interface IMs1FeatureFile
    {
        public IEnumerable<ISingleChargeMs1Feature> GetMs1Features();
    }
}
