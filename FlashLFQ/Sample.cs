using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    /// <summary>
    /// Defines the type of sample this is. 
    /// Types include MS1 label-free, MS1 labeled (e.g., SILAC), MS2+ label-free (e.g., DIA), MS2+ labeled (TMT).
    /// The type defines how this sample will be quantified.
    /// </summary>
    public enum SampleType { Ms1LabelFree, Ms1Labeled, MsnLabelFree, MsnLabeled }

    /// <summary>
    /// Defines information for a sample. A single mass spec data file can have multiple samples (TMT, SILAC),
    /// and a sample can be spread across multiple spectra files (fractionated samples). In unfractionated LFQ,
    /// a sample is synonymous with an MS data file.
    /// </summary>
    public class Sample
    {
        public readonly string SampleGroup;
        public readonly string SampleId;
        public readonly List<string> FractionIds;
        public readonly SampleType SampleType;

        public Sample(string condition, string sampleId, List<string> fractionIds, SampleType sampleType)
        {
            this.SampleGroup = condition;
            this.SampleId = sampleId;
            this.FractionIds = fractionIds;
            this.SampleType = sampleType;
        }

        public override int GetHashCode()
        {
            return (SampleGroup + SampleId).GetHashCode();
        }

        public override bool Equals(object obj)
        {
            var otherSample = (Sample)obj;

            if (otherSample == null)
            {
                return false;
            }

            return otherSample.SampleGroup == this.SampleGroup && otherSample.SampleId == this.SampleId;
        }
    }
}
