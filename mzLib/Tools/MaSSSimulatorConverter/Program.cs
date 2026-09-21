using Readers.MaSSSimulator;

if (args.Length != 3)
{
    Console.Error.WriteLine("Usage: MaSSSimulatorConverter <spectra> <truth.rst> <output.mgf>");
    return 2;
}

var spectra = new MaSSSimulatorSpectrumFile(args[0]);
spectra.ApplyTruth(args[1]);
spectra.WriteMgf(args[2]);
Console.WriteLine($"Converted {spectra.Results.Count} spectra to {args[2]}");
return 0;
