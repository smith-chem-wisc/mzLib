using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public interface ISingleChargeMs1Feature
    {
        double Mz { get; }
        int Charge { get; }
        double RetentionTimeStart { get; }
        double RetentionTimeEnd { get; }
        double Intensity { get; }
        int? NumberOfIsotopes { get; }
        double? FractionalIntensity { get; }
    }
}
