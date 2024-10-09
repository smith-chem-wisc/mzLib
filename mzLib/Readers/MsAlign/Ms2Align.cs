using MassSpectrometry;
namespace Readers;

public class Ms2Align : MsAlign
{
    public override int DefaultMsnOrder => 2;
    public Ms2Align(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile) { }
    public Ms2Align(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile) { }
    public Ms2Align(string filePath) : base(filePath) { }
}