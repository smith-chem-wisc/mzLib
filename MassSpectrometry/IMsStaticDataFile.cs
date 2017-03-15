namespace MassSpectrometry
{
    public interface IMsStaticDataFile<TScan>
        where TScan : IMsDataScan<IMzSpectrum<IMzPeak>>
    {
    }
}