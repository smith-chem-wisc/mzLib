using Chemistry;
using System;
using System.Collections.Generic;
using System.Text;

namespace FlashLFQ
{
    public class SampleComponent
    {
        public SpectraFileInfo SpectraFileInfo;
        public readonly string SampleGroup;
        public readonly int Sample;
        public readonly int Fraction;
        public readonly int Replicate;

        public readonly SampleType SampleType;
        public readonly ChemicalLabel ChemicalLabel;

        public SampleComponent(SpectraFileInfo file, string sampleGroup, int sample, int fraction, int replicate, SampleType sampleType, ChemicalLabel label = null)
        {
            
        }
    }
}
