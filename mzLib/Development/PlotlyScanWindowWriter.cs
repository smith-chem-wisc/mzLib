using System.Globalization;
using System.Text;
using MassSpectrometry;

namespace Development;

internal static class PlotlyScanWindowWriter
{
    internal static string BuildHtml(IReadOnlyList<MsDataScan> scans, string rawPath, double rtCenter, double windowMinutes)
    {
        var xs = new StringBuilder();
        var ys = new StringBuilder();
        var intensities = new StringBuilder();
        var labels = new StringBuilder();

        var first = true;
        foreach (var scan in scans)
        {
            var xArray = scan.MassSpectrum.XArray;
            var yArray = scan.MassSpectrum.YArray;
            for (int i = 0; i < xArray.Length; i++)
            {
                if (!first)
                {
                    xs.Append(',');
                    ys.Append(',');
                    intensities.Append(',');
                    labels.Append(',');
                }

                first = false;
                xs.Append(xArray[i].ToString(CultureInfo.InvariantCulture));
                ys.Append(scan.RetentionTime.ToString(CultureInfo.InvariantCulture));
                intensities.Append(yArray[i].ToString(CultureInfo.InvariantCulture));
                labels.Append('"').Append($"scan {scan.OneBasedScanNumber}").Append('"');
            }
        }

        var title = $"MS1 window around RT {rtCenter.ToString(CultureInfo.InvariantCulture)} min (+/- {windowMinutes.ToString(CultureInfo.InvariantCulture)} min) from {Path.GetFileName(rawPath)}";
        var html = new StringBuilder();
        html.AppendLine("<!doctype html>");
        html.AppendLine("<html>");
        html.AppendLine("<head>");
        html.AppendLine("  <meta charset=\"utf-8\" />");
        html.AppendLine("  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>");
        html.AppendLine("  <style>body { margin: 0; font-family: sans-serif; } #plot { width: 100vw; height: 100vh; }</style>");
        html.AppendLine("</head>");
        html.AppendLine("<body>");
        html.AppendLine("  <div id=\"plot\"></div>");
        html.AppendLine("  <script>");
        html.AppendLine($"    const trace = {{ x: [{xs}], y: [{ys}], mode: 'markers', type: 'scattergl', marker: {{ size: 5, color: [{intensities}], colorscale: 'Viridis', showscale: true }}, text: [{labels}], hovertemplate: '%{{text}}<br>m/z=%{{x}}<br>RT=%{{y}}<br>Intensity=%{{marker.color}}<extra></extra>' }};");
        html.AppendLine($"    Plotly.newPlot('plot', [trace], {{ title: '{title.Replace("'", "\\'")}', xaxis: {{ title: 'm/z' }}, yaxis: {{ title: 'Retention time (min)' }}, margin: {{ l: 60, r: 20, t: 60, b: 60 }} }});");
        html.AppendLine("  </script>");
        html.AppendLine("</body>");
        html.AppendLine("</html>");

        return html.ToString();
    }

    internal static void WriteHtml(IReadOnlyList<MsDataScan> scans, string rawPath, double rtCenter, double windowMinutes, string outputHtml)
    {
        File.WriteAllText(outputHtml, BuildHtml(scans, rawPath, rtCenter, windowMinutes));
    }
}
