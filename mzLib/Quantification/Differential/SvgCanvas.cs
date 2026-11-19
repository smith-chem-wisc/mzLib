using System.Globalization;
using System.Text;

namespace Quantification.Differential;

/// <summary>
/// A minimal SVG drawing surface for the statistics report's figures, written by our own code rather than a plotting
/// library (STATS-FRAMEWORK section 6). Every number goes through <see cref="N"/>: invariant culture, two decimals, no
/// negative zero, so the same figure gives the same bytes on every machine. Lines end in <c>\n</c>.
/// </summary>
internal sealed class SvgCanvas
{
    private readonly StringBuilder _body = new();

    internal SvgCanvas(double width, double height)
    {
        Width = width;
        Height = height;
    }

    internal double Width { get; }
    internal double Height { get; }

    /// <summary>A coordinate or size: invariant culture, at most two decimals, never <c>-0</c>.</summary>
    internal static string N(double value)
    {
        double rounded = Math.Round(value, 2, MidpointRounding.AwayFromZero);
        if (rounded == 0) rounded = 0;
        return rounded.ToString("0.##", CultureInfo.InvariantCulture);
    }

    /// <summary>A tick or cell label with a fixed number of decimals, invariant culture, never <c>-0</c>.</summary>
    internal static string Label(double value, int decimals)
    {
        double rounded = Math.Round(value, decimals, MidpointRounding.AwayFromZero);
        if (rounded == 0) rounded = 0;
        return rounded.ToString("F" + decimals.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
    }

    internal static string Escape(string text) =>
        text.Replace("&", "&amp;").Replace("<", "&lt;").Replace(">", "&gt;").Replace("\"", "&quot;");

    internal void Rect(double x, double y, double width, double height, string fill, string? stroke = null) =>
        _body.Append($"  <rect x=\"{N(x)}\" y=\"{N(y)}\" width=\"{N(Math.Max(0, width))}\" height=\"{N(Math.Max(0, height))}\" " +
                     $"fill=\"{fill}\"{(stroke is null ? "" : $" stroke=\"{stroke}\"")}/>\n");

    internal void Line(double x1, double y1, double x2, double y2, string stroke, double width = 1, bool dashed = false) =>
        _body.Append($"  <line x1=\"{N(x1)}\" y1=\"{N(y1)}\" x2=\"{N(x2)}\" y2=\"{N(y2)}\" stroke=\"{stroke}\" " +
                     $"stroke-width=\"{N(width)}\"{(dashed ? " stroke-dasharray=\"4 3\"" : "")}/>\n");

    internal void Circle(double cx, double cy, double r, string fill, double opacity = 1) =>
        _body.Append($"  <circle cx=\"{N(cx)}\" cy=\"{N(cy)}\" r=\"{N(r)}\" fill=\"{fill}\"" +
                     $"{(opacity < 1 ? $" fill-opacity=\"{N(opacity)}\"" : "")}/>\n");

    internal void Polygon(IEnumerable<(double X, double Y)> points, string fill) =>
        _body.Append($"  <polygon points=\"{string.Join(" ", points.Select(p => $"{N(p.X)},{N(p.Y)}"))}\" fill=\"{fill}\"/>\n");

    internal void Polyline(IEnumerable<(double X, double Y)> points, string stroke, double width = 1.5) =>
        _body.Append($"  <polyline points=\"{string.Join(" ", points.Select(p => $"{N(p.X)},{N(p.Y)}"))}\" fill=\"none\" " +
                     $"stroke=\"{stroke}\" stroke-width=\"{N(width)}\"/>\n");

    internal void Text(double x, double y, string text, string anchor = "start", double size = 12, double rotate = 0)
    {
        string transform = rotate == 0 ? "" : $" transform=\"rotate({N(rotate)} {N(x)} {N(y)})\"";
        string font = size == 12 ? "" : $" font-size=\"{N(size)}\"";
        _body.Append($"  <text x=\"{N(x)}\" y=\"{N(y)}\" text-anchor=\"{anchor}\"{font}{transform}>{Escape(text)}</text>\n");
    }

    internal byte[] ToBytes()
    {
        var svg = new StringBuilder();
        svg.Append($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{N(Width)}\" height=\"{N(Height)}\" " +
                   $"viewBox=\"0 0 {N(Width)} {N(Height)}\" font-family=\"sans-serif\" font-size=\"12\">\n");
        svg.Append("  <rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        svg.Append(_body);
        svg.Append("</svg>\n");
        return new UTF8Encoding(false).GetBytes(svg.ToString());
    }
}

/// <summary>A linear map from data values to SVG coordinates.</summary>
internal readonly record struct SvgScale(double DomainLo, double DomainHi, double RangeLo, double RangeHi)
{
    internal double Map(double value) => RangeLo + (value - DomainLo) / (DomainHi - DomainLo) * (RangeHi - RangeLo);
}

/// <summary>An axis: its domain, its tick values and how many decimals the tick labels carry.</summary>
internal sealed record SvgAxis(double Lo, double Hi, IReadOnlyList<double> Ticks, int Decimals)
{
    /// <summary>
    /// An axis covering [<paramref name="min"/>, <paramref name="max"/>] with about <paramref name="target"/> ticks at a
    /// step of 1, 2 or 5 times a power of ten; the domain is widened to whole steps. A degenerate range is widened.
    /// </summary>
    internal static SvgAxis Nice(double min, double max, int target = 5)
    {
        if (!double.IsFinite(min) || !double.IsFinite(max)) (min, max) = (0, 1);
        if (max <= min)
        {
            double pad = min != 0 ? Math.Abs(min) * 0.1 : 1;
            (min, max) = (min - pad, max + pad);
        }
        double raw = (max - min) / target;
        double magnitude = Math.Pow(10, Math.Floor(Math.Log10(raw)));
        double normalized = raw / magnitude;
        double step = (normalized < 1.5 ? 1 : normalized < 3 ? 2 : normalized < 7 ? 5 : 10) * magnitude;
        double lo = Math.Floor(min / step) * step, hi = Math.Ceiling(max / step) * step;
        int count = (int)Math.Round((hi - lo) / step);
        var ticks = Enumerable.Range(0, count + 1).Select(k => lo + k * step).ToArray();
        int decimals = Math.Max(0, -(int)Math.Floor(Math.Log10(step) + 1e-9));
        return new SvgAxis(lo, hi, ticks, decimals);
    }
}

/// <summary>The shared plot frame: title, axes, ticks and axis labels.</summary>
internal static class SvgPlot
{
    internal const double Left = 72, Right = 20, Top = 40, Bottom = 56;
    internal const string Ink = "#333333";
    internal const string Points = "#2166ac";
    internal const string Accent = "#d6604d";
    internal const string Bars = "#4c78a8";

    internal static (SvgScale X, SvgScale Y) Frame(SvgCanvas c, string title, string xLabel, string yLabel, SvgAxis x, SvgAxis y)
    {
        var sx = new SvgScale(x.Lo, x.Hi, Left, c.Width - Right);
        var sy = new SvgScale(y.Lo, y.Hi, c.Height - Bottom, Top);
        c.Text(c.Width / 2, 22, title, "middle", 14);
        c.Line(Left, c.Height - Bottom, c.Width - Right, c.Height - Bottom, Ink);
        c.Line(Left, Top, Left, c.Height - Bottom, Ink);
        foreach (double t in x.Ticks)
        {
            double px = sx.Map(t);
            c.Line(px, c.Height - Bottom, px, c.Height - Bottom + 5, Ink);
            c.Text(px, c.Height - Bottom + 18, SvgCanvas.Label(t, x.Decimals), "middle");
        }
        foreach (double t in y.Ticks)
        {
            double py = sy.Map(t);
            c.Line(Left - 5, py, Left, py, Ink);
            c.Text(Left - 8, py + 4, SvgCanvas.Label(t, y.Decimals), "end");
        }
        c.Text((Left + c.Width - Right) / 2, c.Height - 14, xLabel, "middle");
        c.Text(18, (Top + c.Height - Bottom) / 2, yLabel, "middle", rotate: -90);
        return (sx, sy);
    }
}
