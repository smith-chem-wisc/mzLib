using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class SampleComponent
    {
        public readonly string SampleGroup;
        public readonly int Sample;
        public readonly int Fraction;
        public readonly int Replicate;

        public readonly SampleType SampleType;
        public string Label;
        public double LabelMass;

        public SampleComponent()
        {

        }
    }
}
