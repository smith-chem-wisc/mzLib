using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Omics.BioPolymerGroup;

namespace Quantification.Differential;

/// <summary>A pair of samples' Pearson correlation over the features both have: <see cref="R"/> is NaN below 3.</summary>
internal sealed record SampleCorrelation(string A, string B, double R, int N);

/// <summary>One feature on limma's plotSA scale: average log2 intensity, sigma^(1/2) and the prior's (s0^2)^(1/4).</summary>
internal sealed record MeanVariancePoint(string FeatureId, double X, double Y, double PriorY, double ResidualVariance, double PriorVariance);

/// <summary>
/// The statistics report's diagnostic figures (STATS-FRAMEWORK section 6, item 4): the numbers behind each, and each
/// drawn as an SVG with the table of its data. Each table is written as the results table is: a <see cref="TsvColumn{T}"/>
/// schema through <see cref="TsvWriter"/>, numbers by <see cref="DifferentialColumns.Number"/> and text by
/// <see cref="DifferentialResultWriter.CleanText"/>, so <c>DEF-DIFF-FILE</c>'s encoding holds: tab-separated, <c>\n</c>
/// line ends, round-trip numbers in invariant culture, a missing value as an empty cell, and a tab, CR or LF inside text
/// written as a space. The metadata's <c>cleaned_text_cells</c> counts the results table's cells only; a figure-table cell
/// cleaned by the same rule is not counted there.
/// </summary>
internal static class DiagnosticFigures
{
    internal const int PValueBinCount = 20;

    private static readonly TsvColumn<(double Low, double High, int Features)>[] PValueTable =
    {
        new("bin_low", b => DifferentialColumns.Number(b.Low)!),
        new("bin_high", b => DifferentialColumns.Number(b.High)!),
        new("features", b => DifferentialColumns.Integer(b.Features)!),
    };

    private static readonly TsvColumn<(string FeatureId, double Log2, double P)>[] VolcanoTable =
    {
        new("feature_id", p => DifferentialResultWriter.CleanText(p.FeatureId)!),
        new("log2_effect", p => DifferentialColumns.Number(p.Log2)!),
        new("p_value", p => DifferentialColumns.Number(p.P)!),
        new("neg_log10_p", p => DifferentialColumns.Number(p.P > 0 ? -Math.Log10(p.P) : null)!),
    };

    private static readonly TsvColumn<MeanVariancePoint>[] MeanVarianceTable =
    {
        new("feature_id", p => DifferentialResultWriter.CleanText(p.FeatureId)!),
        new("average_log2", p => DifferentialColumns.Number(p.X)!),
        new("residual_variance", p => DifferentialColumns.Number(p.ResidualVariance)!),
        new("sqrt_sigma", p => DifferentialColumns.Number(p.Y)!),
        new("prior_variance", p => DifferentialColumns.Number(p.PriorVariance)!),
        new("prior_sqrt_sigma", p => DifferentialColumns.Number(p.PriorY)!),
    };

    private static readonly TsvColumn<SampleCorrelation>[] CorrelationTable =
    {
        new("sample_a", p => DifferentialResultWriter.CleanText(p.A)!),
        new("sample_b", p => DifferentialResultWriter.CleanText(p.B)!),
        new("r", p => DifferentialColumns.Number(p.R)!),
        new("features", p => DifferentialColumns.Integer(p.N)!),
    };

    private static readonly TsvColumn<(string Contrast, DifferentialGlobalShift Shift)>[] GlobalShiftTable =
    {
        new("contrast_id", e => DifferentialResultWriter.CleanText(e.Contrast)!),
        new("stratum", e => DifferentialResultWriter.CleanText(e.Shift.Stratum)!),
        new("quant_basis", e => DifferentialResultWriter.CleanText(e.Shift.QuantBasis)!),
        new("value", e => DifferentialColumns.Number(e.Shift.Value)!),
        new("peptides", e => DifferentialColumns.Integer(e.Shift.Peptides)!),
    };

    /// <summary>
    /// Counts p-values in 20 bins of width 0.05. Bin k holds 0.05k &lt;= p &lt; 0.05(k + 1), computed as floor(20p); p = 1
    /// goes in the last bin. Non-finite values are not counted.
    /// </summary>
    internal static int[] PValueBins(IEnumerable<double> pValues)
    {
        var bins = new int[PValueBinCount];
        foreach (double p in pValues.Where(double.IsFinite))
            bins[Math.Clamp((int)Math.Floor(p * PValueBinCount), 0, PValueBinCount - 1)]++;
        return bins;
    }

    /// <summary>
    /// Pearson correlation of each pair of samples (in sample order, first index lower), over the features where both
    /// samples have a finite value. With fewer than 3 such features, r is NaN.
    /// </summary>
    internal static IReadOnlyList<SampleCorrelation> Correlations(SampleLog2 samples)
    {
        var result = new List<SampleCorrelation>();
        for (int a = 0; a < samples.SampleIds.Count; a++)
            for (int b = a + 1; b < samples.SampleIds.Count; b++)
            {
                var pairs = samples.Values[a].Zip(samples.Values[b])
                    .Where(p => double.IsFinite(p.First) && double.IsFinite(p.Second)).ToList();
                double r = double.NaN;
                if (pairs.Count >= 3)
                {
                    double mx = pairs.Average(p => p.First), my = pairs.Average(p => p.Second);
                    double sxy = 0, sxx = 0, syy = 0;
                    foreach (var (x, y) in pairs)
                    {
                        sxy += (x - mx) * (y - my);
                        sxx += (x - mx) * (x - mx);
                        syy += (y - my) * (y - my);
                    }
                    r = sxx > 0 && syy > 0 ? Math.Clamp(sxy / Math.Sqrt(sxx * syy), -1, 1) : double.NaN;
                }
                result.Add(new SampleCorrelation(samples.SampleIds[a], samples.SampleIds[b], r, pairs.Count));
            }
        return result;
    }

    /// <summary>
    /// limma's <c>plotSA</c> scale: x = average log2 intensity, y = sigma^(1/2) = (residual variance)^(1/4), and the prior
    /// drawn at (prior variance)^(1/4). Features with a non-finite x or y are left out.
    /// </summary>
    internal static IReadOnlyList<MeanVariancePoint> MeanVariancePoints(DiagnosticFit fit) =>
        Enumerable.Range(0, fit.FeatureIds.Count)
            .Select(i => new MeanVariancePoint(fit.FeatureIds[i], fit.AverageLog2[i], Math.Pow(fit.ResidualVariance[i], 0.25),
                Math.Pow(fit.PriorVariance[i], 0.25), fit.ResidualVariance[i], fit.PriorVariance[i]))
            .Where(p => double.IsFinite(p.X) && double.IsFinite(p.Y))
            .ToList();

    internal static (byte[] Svg, byte[] Tsv) PValueHistogram(string title, IReadOnlyList<double> pValues)
    {
        var bins = PValueBins(pValues);
        var c = new SvgCanvas(640, 480);
        var x = new SvgAxis(0, 1, new[] { 0, 0.2, 0.4, 0.6, 0.8, 1.0 }, 1);
        var (sx, sy) = SvgPlot.Frame(c, title, "p-value", "features", x, SvgAxis.Nice(0, Math.Max(1, bins.Max())));
        for (int k = 0; k < PValueBinCount; k++)
        {
            double x0 = sx.Map(k / (double)PValueBinCount), x1 = sx.Map((k + 1) / (double)PValueBinCount);
            c.Rect(x0 + 0.5, sy.Map(bins[k]), x1 - x0 - 1, sy.Map(0) - sy.Map(bins[k]), SvgPlot.Bars);
        }
        return (c.ToBytes(), Table(PValueTable, Enumerable.Range(0, PValueBinCount)
            .Select(k => (k / (double)PValueBinCount, (k + 1) / (double)PValueBinCount, bins[k]))));
    }

    internal static (byte[] Svg, byte[] Tsv) Volcano(string title, IReadOnlyList<(string FeatureId, double Log2, double P)> points)
    {
        var drawn = points.Where(p => double.IsFinite(p.Log2) && double.IsFinite(p.P)).ToList();
        double spread = Math.Max(1e-9, drawn.Count == 0 ? 1 : drawn.Max(p => Math.Abs(p.Log2)));
        double top = drawn.Where(p => p.P > 0).Select(p => -Math.Log10(p.P)).DefaultIfEmpty(1).Max();
        var c = new SvgCanvas(640, 480);
        var (sx, sy) = SvgPlot.Frame(c, title, "log2 effect", "-log10(p)", SvgAxis.Nice(-spread, spread), SvgAxis.Nice(0, Math.Max(1, top)));
        c.Line(sx.Map(0), SvgPlot.Top, sx.Map(0), c.Height - SvgPlot.Bottom, "#bbbbbb", dashed: true);
        foreach (var p in drawn)
        {
            if (p.P > 0)
            {
                c.Circle(sx.Map(p.Log2), sy.Map(-Math.Log10(p.P)), 3, SvgPlot.Points, 0.6);
                continue;
            }
            double px = sx.Map(p.Log2);
            c.Polygon(new[] { (px, SvgPlot.Top), (px - 4, SvgPlot.Top + 7), (px + 4, SvgPlot.Top + 7) }, SvgPlot.Accent);
        }
        return (c.ToBytes(), Table(VolcanoTable, points));
    }

    internal static (byte[] Svg, byte[] Tsv) MeanVariance(string title, DiagnosticFit fit)
    {
        var points = MeanVariancePoints(fit);
        var ys = points.Select(p => p.Y).Concat(points.Select(p => p.PriorY).Where(double.IsFinite)).ToList();
        var c = new SvgCanvas(640, 480);
        var (sx, sy) = SvgPlot.Frame(c, title, "average log2 intensity", "sqrt(sigma)",
            SvgAxis.Nice(points.Select(p => p.X).DefaultIfEmpty(0).Min(), points.Select(p => p.X).DefaultIfEmpty(1).Max()),
            SvgAxis.Nice(ys.DefaultIfEmpty(0).Min(), ys.DefaultIfEmpty(1).Max()));
        foreach (var p in points)
            c.Circle(sx.Map(p.X), sy.Map(p.Y), 3, SvgPlot.Points, 0.6);
        var prior = points.Where(p => double.IsFinite(p.PriorY)).OrderBy(p => p.X).ThenBy(p => p.FeatureId, StringComparer.Ordinal)
            .Select(p => (sx.Map(p.X), sy.Map(p.PriorY))).ToList();
        if (prior.Count > 1) c.Polyline(prior, SvgPlot.Accent, 2);
        return (c.ToBytes(), Table(MeanVarianceTable, points));
    }

    internal static (byte[] Svg, byte[] Tsv) SampleCorrelationHeatmap(string title, SampleLog2 samples)
    {
        var pairs = Correlations(samples);
        int n = samples.SampleIds.Count;
        const double cell = 28, left = 110, top = 110;
        double width = Math.Max(420, left + n * cell + 30), height = top + n * cell + 90;
        var c = new SvgCanvas(width, height);
        c.Text(width / 2, 22, title, "middle", 14);
        double Lookup(int a, int b)
        {
            if (a == b) return 1;
            var (i, j) = a < b ? (a, b) : (b, a);
            return pairs.Single(p => p.A == samples.SampleIds[i] && p.B == samples.SampleIds[j]).R;
        }
        for (int a = 0; a < n; a++)
        {
            c.Text(left - 6, top + a * cell + cell / 2 + 4, samples.SampleIds[a], "end");
            c.Text(left + a * cell + cell / 2, top - 6, samples.SampleIds[a], "start", rotate: -45);
            for (int b = 0; b < n; b++)
            {
                double r = Lookup(a, b);
                c.Rect(left + b * cell, top + a * cell, cell, cell, Color(r), "#ffffff");
                if (n <= 16 && double.IsFinite(r))
                    c.Text(left + b * cell + cell / 2, top + a * cell + cell / 2 + 4, SvgCanvas.Label(r, 2), "middle", 9);
            }
        }
        double legendY = top + n * cell + 30;
        c.Text(left, legendY - 8, "Pearson r (white: fewer than 3 shared features)", "start", 11);
        double[] stops = { -1, -0.5, 0, 0.5, 1 };
        for (int k = 0; k < stops.Length; k++)
        {
            c.Rect(left + k * 44, legendY, 40, 14, Color(stops[k]), "#999999");
            c.Text(left + k * 44 + 20, legendY + 28, SvgCanvas.Label(stops[k], 1), "middle", 11);
        }
        return (c.ToBytes(), Table(CorrelationTable, pairs));
    }

    internal static (byte[] Svg, byte[] Tsv) GlobalShift(string title, IReadOnlyList<(string Contrast, DifferentialGlobalShift Shift)> entries)
    {
        const double left = 240, right = 150, top = 56, rowHeight = 28;
        double width = 720, height = top + entries.Count * rowHeight + 60;
        var values = entries.Select(e => e.Shift.Value).Where(double.IsFinite).ToList();
        double spread = Math.Max(0.1, values.Select(v => Math.Abs(v)).DefaultIfEmpty(0.1).Max());
        var axis = SvgAxis.Nice(-spread, spread);
        var sx = new SvgScale(axis.Lo, axis.Hi, left, width - right);
        var c = new SvgCanvas(width, height);
        c.Text(width / 2, 22, title, "middle", 14);
        double bottom = top + entries.Count * rowHeight;
        c.Line(sx.Map(0), top - 8, sx.Map(0), bottom, "#bbbbbb", dashed: true);
        c.Line(left, bottom, width - right, bottom, SvgPlot.Ink);
        foreach (double t in axis.Ticks)
        {
            c.Line(sx.Map(t), bottom, sx.Map(t), bottom + 5, SvgPlot.Ink);
            c.Text(sx.Map(t), bottom + 18, SvgCanvas.Label(t, axis.Decimals), "middle");
        }
        c.Text((left + width - right) / 2, bottom + 40, "global shift (log2, numerator over denominator, before normalization)", "middle");
        for (int k = 0; k < entries.Count; k++)
        {
            var (contrast, shift) = entries[k];
            double y = top + k * rowHeight + rowHeight / 2;
            c.Text(left - 10, y + 4, $"{contrast}, {shift.Stratum}{(shift.QuantBasis is null ? "" : ", " + shift.QuantBasis)}", "end");
            string peptides = $"{shift.Peptides.ToString(CultureInfo.InvariantCulture)} peptides";
            if (double.IsFinite(shift.Value))
                c.Circle(sx.Map(shift.Value), y, 5, SvgPlot.Points);
            else
                peptides = "no value (" + peptides + ")";
            c.Text(width - right + 10, y + 4, peptides, "start", 11);
        }
        return (c.ToBytes(), Table(GlobalShiftTable, entries));
    }

    /// <summary>Red through white to blue for r in [-1, 1]; white for a missing r.</summary>
    internal static string Color(double r)
    {
        if (!double.IsFinite(r)) return "#ffffff";
        (int R, int G, int B) neg = (0xb2, 0x18, 0x2b), mid = (0xf7, 0xf7, 0xf7), pos = (0x21, 0x66, 0xac);
        var (from, to, t) = r < 0 ? (mid, neg, -r) : (mid, pos, r);
        t = Math.Clamp(t, 0, 1);
        int Mix(int a, int b) => (int)Math.Round(a + (b - a) * t, MidpointRounding.AwayFromZero);
        return $"#{Mix(from.R, to.R):x2}{Mix(from.G, to.G):x2}{Mix(from.B, to.B):x2}";
    }

    /// <summary>
    /// A figure's table, written as <see cref="DifferentialResultWriter"/> writes the results table: the header and each row
    /// from <see cref="TsvWriter"/>, each line ended by <c>\n</c>, as UTF-8 without a byte-order mark.
    /// </summary>
    private static byte[] Table<T>(IReadOnlyList<TsvColumn<T>> schema, IEnumerable<T> rows)
    {
        var text = new StringBuilder(TsvWriter.HeaderLine(schema)).Append('\n');
        foreach (var row in rows)
            text.Append(TsvWriter.RowLine(schema, row)).Append('\n');
        return new UTF8Encoding(false).GetBytes(text.ToString());
    }
}
