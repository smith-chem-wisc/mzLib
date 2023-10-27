using MassSpectrometry;

namespace Readers;

public class GenericMsDataFile : MsDataFile
{
    public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
    {
        throw new NotImplementedException();
    }
    public override SourceFile GetSourceFile()
    {
        throw new NotImplementedException();
    }

    public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
    {
        throw new NotImplementedException();
    }

    public override void CloseDynamicConnection()
    {
        throw new NotImplementedException();
    }

    public override void InitiateDynamicConnection()
    {
        throw new NotImplementedException();
    }
    public GenericMsDataFile(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
    {

    }
    public GenericMsDataFile(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
    {

    }

    public GenericMsDataFile(string filePath) : base(filePath)
    {

    }
}