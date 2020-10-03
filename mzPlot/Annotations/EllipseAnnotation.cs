using System;
using System.Collections.Generic;
using System.Text;
using OxyPlot;

namespace mzPlot
{
    public class EllipseAnnotation : mzPlotAnnotation
    {
        public EllipseAnnotation() : base(AnnotationTypes.Data, 0, 0)
        {

        }

        public override void Render(IRenderContext rc)
        {
            throw new NotImplementedException();
        }
    }
}
