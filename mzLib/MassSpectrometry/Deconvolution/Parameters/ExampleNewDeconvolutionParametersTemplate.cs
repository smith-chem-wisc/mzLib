﻿using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    public class ExampleNewDeconvolutionParametersTemplate : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.ExampleNewDeconvolutionTemplate;
        public ExampleNewDeconvolutionParametersTemplate(int minCharge, int maxCharge, Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {

        }
    }
}
