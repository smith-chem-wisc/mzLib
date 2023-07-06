namespace Proteomics.ProteolyticDigestion
{
    public enum CleavageSpecificity
    {
        None,
        Semi,
        Full,
        SingleN, // 5' for RNA
        SingleC, // 3' for 
        Unknown //used for fast Semi/NonSpecific searches when peptide is cleaved post-search
    }
}