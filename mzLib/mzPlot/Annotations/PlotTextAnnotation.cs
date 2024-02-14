using OxyPlot;
using OxyPlot.Annotations;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class PlotTextAnnotation : Annotation
    {
        /// <summary>
        /// Creates a "floating" text annotation, i.e., the annotation's x,y position is coupled to the chart area and not a data point.
        /// </summary>
        public PlotTextAnnotation()
        {

        }

        public string Text { get; set; }

        /// <summary>
        /// The x-position of the annotation, in pixels. Relative to the left of the chart (higher X is more to the right). Can be negative.
        /// </summary>
        public double X { get; set; }

        /// <summary>
        /// The y-position of the annotation, in pixels. Relative to the top of the chart (higher Y is lower on the chart). Can be negative.
        /// </summary>
        public double Y { get; set; }

        public override void Render(IRenderContext rc)
        {
            base.Render(rc);
            double pX = PlotModel.PlotArea.Left + X;
            double pY = PlotModel.PlotArea.Top + Y;
            rc.DrawMultilineText(new ScreenPoint(pX, pY), Text, TextColor, Font, FontSize, FontWeight);
        }
    }
}
