using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil
{
    /// <summary>
    /// PpmToleranceWithNotch is only used for neutral mass, considering possible deconvolution errors
    /// </summary>
    public class PpmToleranceWithNotch : Tolerance
    {
        private readonly double[] AcceptableSortedMassShifts;
        public const double NotchStep = 1.00335483810; //C13MinusC12 = 1.00335483810
        public PpmTolerance PpmTolerance { get; private set; }

        public PpmToleranceWithNotch(double value, int positiveNotch, int negativeNotch)
            : base(value)
        {
            var massShifts = new List<double> { 0 };
            for (int i = 1; i <= positiveNotch; i++)
            {
                massShifts.Add(NotchStep * i);
            }
            for (int i = 1; i <= negativeNotch; i++)
            {
                massShifts.Add(-NotchStep * i);
            }
            AcceptableSortedMassShifts = massShifts.OrderBy(p => p).ToArray();
            PpmTolerance = new PpmTolerance(value);
        }

        public override double GetMaximumValue(double mean)
        {
            return (PpmTolerance.GetMaximumValue(AcceptableSortedMassShifts[AcceptableSortedMassShifts.Length - 1] + mean));
        }

        public override double GetMinimumValue(double mean)
        {
            return (PpmTolerance.GetMinimumValue(AcceptableSortedMassShifts[0] + mean));
        }

        public override DoubleRange GetRange(double mean)
        {
            return new DoubleRange(GetMinimumValue(mean), GetMaximumValue(mean));
        }

        public override bool Within(double experimental, double theoretical)
        {
            for (int i = 0; i < AcceptableSortedMassShifts.Length; i++)
            {
                if (PpmTolerance.Within(experimental, theoretical + AcceptableSortedMassShifts[i]))
                {
                    return true;
                }
            }
            return false;
        }

        public override Tolerance UpdateTolerance(double newValue)
        {
            this.PpmTolerance = new PpmTolerance(newValue);
            return this;
        }
    }
}
