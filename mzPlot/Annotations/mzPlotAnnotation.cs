using OxyPlot;
using OxyPlot.Annotations;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public abstract class mzPlotAnnotation : Annotation
    {
        public enum AnnotationTypes { Floating, Data }

        public readonly AnnotationTypes AnnotationType;
        public readonly double X;
        public readonly double Y;

        public mzPlotAnnotation(AnnotationTypes annotationType, double x, double y)
        {
            this.AnnotationType = annotationType;
            this.X = x;
            this.Y = y;
        }
    }
}
