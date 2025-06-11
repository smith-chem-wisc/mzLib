using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil
{
    public class PpmToleranceWithNotch : Tolerance
    {
        private readonly double[] AcceptableSortedMassShifts;
        public const double NotchStep = 1.00335483810; //C13MinusC12 = 1.00335483810
        public PpmTolerance PpmTolerance { get; init; }

        public PpmToleranceWithNotch(double value, int protonNotches)
            : base(value)
        {
            var massShifts = new List<double> { 0 };
            for (int i = 1; i <= protonNotches; i++)
            {
                massShifts.Add(NotchStep * i);
                massShifts.Add(-NotchStep * i);
            }
            AcceptableSortedMassShifts = massShifts.OrderBy(Math.Abs).ThenBy(p => p).ToArray();
            PpmTolerance = new PpmTolerance(value);
        }

        public PpmToleranceWithNotch(double value, List<int> massShifts)
            : base(value)
        {
            var shifts = new List<double> { 0 };
            foreach (int shift in massShifts)
            {
                shifts.Add(NotchStep * shift);
            }
            AcceptableSortedMassShifts = shifts.OrderBy(Math.Abs).ThenBy(p => p).Distinct().ToArray();
            PpmTolerance = new PpmTolerance(value);
        }

        public override double GetMaximumValue(double mean)
        {
            return (PpmTolerance.GetMaximumValue(AcceptableSortedMassShifts[AcceptableSortedMassShifts.Length - 1] + mean));
        }

        public override double GetMinimumValue(double mean)
        {
            return (PpmTolerance.GetMinimumValue(AcceptableSortedMassShifts[AcceptableSortedMassShifts.Length - 2] + mean));
        }

        public override DoubleRange GetRange(double mean)
        {
            return new DoubleRange(GetMinimumValue(mean), GetMaximumValue(mean));
        }

        public override bool Within(double experimental, double theoretical)
        {
            for (int i = 0; i < AcceptableSortedMassShifts.Length - 1; i++)
            {
                if (PpmTolerance.Within(experimental, theoretical + AcceptableSortedMassShifts[i]))
                {
                    return true;
                }
            }
            return false;
        }
    }
}
