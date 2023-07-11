using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.PeakPicking
{
    public static class PeakPickingExtensions
    {
        public static int GetClosestIndex(this Extremum[] extremaArray, Extremum target,
            ArraySearchOption searchType = ArraySearchOption.Closest)
        {
            int defaultImplementationIndex = Array.BinarySearch(extremaArray, target);
            // Positive index == exact match
            if (defaultImplementationIndex >= 0)
            {
                return defaultImplementationIndex;
            }

            defaultImplementationIndex = ~defaultImplementationIndex; // point to the first value larger than the target.
            if (defaultImplementationIndex == 0)
                return 0;
            if (defaultImplementationIndex == extremaArray.Length)
                defaultImplementationIndex--;

            int closestIndex;
            switch (searchType)
            {
                case ArraySearchOption.Previous:
                    closestIndex = defaultImplementationIndex - 1;
                    break;
                case ArraySearchOption.Closest:
                    double backwardsDiff = target - extremaArray[defaultImplementationIndex - 1];
                    double forwardsDiff = extremaArray[defaultImplementationIndex] - target;
                    closestIndex = backwardsDiff < forwardsDiff
                        ? defaultImplementationIndex - 1
                        : defaultImplementationIndex;
                    break;
                case ArraySearchOption.Next: // Next falls through to default
                default:
                    closestIndex = defaultImplementationIndex;
                    break;
            }

            return closestIndex;
        }

        public static Extremum GetClosestValue(this Extremum[] extremaArray, Extremum target,
            ArraySearchOption searchType = ArraySearchOption.Closest)
        {
            return extremaArray[extremaArray.GetClosestIndex(target, searchType)];
        }
    }
}
