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
        public const double NotchStep = 1.00335483810;
        public PpmTolerance PpmTolerance { get; init; }

        public PpmToleranceWithNotch(double value, int protonNotches)
            : base(value)
        {
            var massShifts = new List<double>();
            for (int i = 1; i <= protonNotches; i++)
            {
                massShifts.Add(NotchStep * i);
                massShifts.Add(-NotchStep * i);
            }
            AcceptableSortedMassShifts = massShifts.OrderBy(p => p).ToArray();
            NumNotches = protonNotches;
            PpmTolerance = new PpmTolerance(value);
        }

        int NumNotches { get; init; }

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
            if (PpmTolerance.Within(experimental, theoretical))
            return true;
            else
            {
                for (int i = 0; i < AcceptableSortedMassShifts.Length - 1; i++)
                {
                    if (PpmTolerance.Within(experimental, theoretical + AcceptableSortedMassShifts[i]))
                    {
                        return true;
                    }
                }
            }
            return false;
        }
    }
}
