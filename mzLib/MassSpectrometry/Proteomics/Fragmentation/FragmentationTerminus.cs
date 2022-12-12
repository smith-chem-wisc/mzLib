namespace Proteomics.Fragmentation
{
    public enum FragmentationTerminus
    {
        Both, //N- and C-terminus
        N, //N-terminus only
        C, //C-terminus only
        None //used for internal fragments, could be used for top down intact mass?
    }
}