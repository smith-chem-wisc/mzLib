using System;
using System.Collections.Generic;
using System.Text;
using OxyPlot;

namespace mzPlot
{
    public class SquareAnnotation : mzPlotAnnotation
    {
        public override void Render(IRenderContext rc)
        {
            throw new NotImplementedException();
        }

        public SquareAnnotation() : base(AnnotationTypes.Data, 0, 0)
        {

        }
    }
}
