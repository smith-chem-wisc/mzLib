using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    public interface IDataReader
    {
        MsDataScan[] Scans { get; }
        SourceFile SourceFile { get; }
        public List<MsDataScan> LoadAllScansFromFile(string filepath, 
            IFilteringParams filteringParams, 
            int maxThreads = -1); 
        public SourceFile GetSourceFile(string filepath);
    }
}