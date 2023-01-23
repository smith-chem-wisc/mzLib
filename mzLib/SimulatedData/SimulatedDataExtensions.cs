using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Plotly.NET.CSharp;
// don't add Plotly.NET as a using. Explicitly refer to it to avoid 
// namespace conflicts.


namespace SimulatedData
{
    public static class SimulatedDataExtensions
    {
        public static Plotly.NET.GenericChart.GenericChart Plot(this SimulatedData simData)
        {
            return Chart.Line<double, double, string>(x: simData.Xarray, y: simData.Yarray)
                .WithXAxisStyle<double, double, string>("m/z")
                .WithYAxisStyle<double, double, string>("Intensity"); 
        }
    }
}
