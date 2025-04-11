using Easy.Common.Extensions;
using Microsoft.ML.Data;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Text;

namespace FlashLFQ.PEP
{
    public class ChromatographicPeakData
    {
        public static readonly IImmutableDictionary<string, string[]> trainingInfos = new Dictionary<string, string[]>
        {
            { "standard", new []
                {
                    "PpmErrorScore",
                    "IntensityScore",
                    "RtScore",
                    "ScanCountScore",
                    "IsotopicDistributionScore",
                    "PpmErrorRaw",
                    "IntensityRaw",
                    "RtPredictionErrorRaw",
                    "ScanCountRaw",
                    "IsotopicPearsonCorrelation"
                } 
            },
            { "reduced", new []
                {
                    "PpmErrorRaw",
                    "IntensityRaw",
                    "RtPredictionErrorRaw",
                    "ScanCountRaw",
                    "IsotopicPearsonCorrelation"
                }
            },
        }.ToImmutableDictionary();

        /// <summary>
        /// These are used for percolator. Trainer must be told the assumed direction for each attribute as it relates to being a true positive
        /// Here, a weight of 1 indicates that the probability of being true is for higher numbers in the set.
        /// A weight of -1 indicates that the probability of being true is for the lower numbers in the set.
        /// </summary>
        public static readonly IImmutableDictionary<string, int> assumedAttributeDirection = new Dictionary<string, int> {
            { "PpmErrorScore", 1 },
            { "IntensityScore", 1 },
            { "RtScore", 1 },
            { "ScanCountScore", 1 },
            { "IsotopicDistributionScore", 1 },
            { "PpmErrorRaw", -1 },
            { "IntensityRaw", 1 },
            { "RtPredictionErrorRaw", -1 },
            { "ScanCountRaw", -1 },
            { "IsotopicPearsonCorrelation", 1 }
            }.ToImmutableDictionary();

        public string ToString(string searchType)
        {
            StringBuilder sb = new StringBuilder();
            var variablesToOutput = ChromatographicPeakData.trainingInfos[searchType];

            foreach (var variable in variablesToOutput)
            {
                var property = typeof(ChromatographicPeakData).GetProperty(variable).GetValue(this, null);
                var floatValue = (float)property;
                sb.Append("\t");
                sb.Append(floatValue.ToString());
            }

            return sb.ToString();
        }

        [LoadColumn(0)]
        public float PpmErrorScore { get; set; }

        [LoadColumn(1)]
        public float IntensityScore { get; set; }

        [LoadColumn(2)]
        public float RtScore { get; set; }

        [LoadColumn(3)]
        public float ScanCountScore { get; set; }

        [LoadColumn(4)]
        public float IsotopicDistributionScore { get; set; }

        [LoadColumn(5)]
        public float PpmErrorRaw { get; set; }

        [LoadColumn(6)]
        public float IntensityRaw { get; set; }

        [LoadColumn(7)]
        public float RtPredictionErrorRaw { get; set; }

        [LoadColumn(8)]
        public float ScanCountRaw { get; set; }

        [LoadColumn(9)]
        public float IsotopicPearsonCorrelation { get; set; }

        [LoadColumn(10)]
        public bool Label { get; set; }

    }
}