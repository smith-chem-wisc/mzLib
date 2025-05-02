using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace FlashLFQ
{
    public enum DonorCriterion
    {
        Score,
        Intensity,
        Neighbors
    }

    public static class DonorCriteriaExtensions
    {
        public static DonorCriterion ParseDonorCriterion(this string donorCriterionString)
        {
            switch(donorCriterionString)
            {
                case ("S"):
                    return DonorCriterion.Score;
                case ("I"):
                    return DonorCriterion.Intensity;
                case ("N"):
                    return DonorCriterion.Neighbors;
                default:
                    throw new MzLibException("Problem with FlashLFQ Parameters: Donor Criterion not recognized");
            }
        }
    }


}
