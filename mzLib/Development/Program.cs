using System.Globalization;
using MassSpectrometry;
using Readers;
using Development;

if (args.Length < 3)
{
    Console.WriteLine("Usage: Development <raw-file> <rt-center-minutes> <output-html> [window-minutes]");
    return;
}

var rawPath = args[0];
var rtCenter = double.Parse(args[1], CultureInfo.InvariantCulture);
var outputHtml = args[2];
var windowMinutes = args.Length > 3 ? double.Parse(args[3], CultureInfo.InvariantCulture) : 1.0;

var reader = MsDataFileReader.GetDataFile(rawPath);
reader.LoadAllStaticData();

var scans = reader.GetAllScansList()
    .Where(s => s.MsnOrder == 1 && s.RetentionTime >= rtCenter - windowMinutes && s.RetentionTime <= rtCenter + windowMinutes)
    .OrderBy(s => s.RetentionTime)
    .ToList();

if (scans.Count == 0)
{
    Console.WriteLine("No MS1 scans found in the requested window.");
    return;
}

PlotlyScanWindowWriter.WriteHtml(scans, rawPath, rtCenter, windowMinutes, outputHtml);
Console.WriteLine($"Wrote {outputHtml}");
