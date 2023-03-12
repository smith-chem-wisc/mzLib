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
        throw new NotSupportedException();
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
    public GenericMsDataFile(int numSpectra) : base(numSpectra, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
    {

    }
    public GenericMsDataFile(MsDataScan[] scans) : base(scans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
    {

    }

    public GenericMsDataFile() : base()
    {

    }

    public GenericMsDataFile(string filePath) : base(filePath)
    {

    }
}