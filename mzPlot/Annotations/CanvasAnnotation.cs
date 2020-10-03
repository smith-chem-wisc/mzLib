using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Controls;

namespace mzPlot
{
    public class CanvasAnnotation : mzPlotAnnotation
    {
        public CanvasAnnotation(Canvas canvas, AnnotationTypes annotationType, double x, double y) : base(annotationType, x, y)
        {
            //canvas.Measure(new Size((int)canvas.Width, 600));
            //canvas.Arrange(new Rect(new Size((int)canvas.Width, 600)));

            //RenderTargetBitmap renderBitmap = new RenderTargetBitmap((int)(canvas.Width), 600, 96, 96, PixelFormats.Pbgra32);

            //renderBitmap.Render(canvas);
            //PngBitmapEncoder encoder = new PngBitmapEncoder();
            //encoder.Frames.Add(BitmapFrame.Create(renderBitmap));
        }
    }
}
