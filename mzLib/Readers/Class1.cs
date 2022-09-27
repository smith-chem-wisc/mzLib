using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    public abstract class DataReader
    {
        public abstract List<MsDataScan> LoadAllScansFromFile(string filepath, 
            IFilteringParams filteringParams, 
            int maxThreads = -1); 
        public abstract SourceFile GetSourceFile(string filepath); 
    }
}