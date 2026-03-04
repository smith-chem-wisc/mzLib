using MassSpectrometry;
namespace Readers;
public class Ms1Align : MsAlign
{
    public override int DefaultMsnOrder => 1;
    public Ms1Align(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile) { }
    public Ms1Align(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile) { }
    public Ms1Align(string filePath) : base(filePath) { }
}