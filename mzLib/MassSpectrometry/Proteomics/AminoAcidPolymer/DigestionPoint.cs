namespace Proteomics.AminoAcidPolymer
{
    public class DigestionPointAndLength
    {
        public DigestionPointAndLength(int index, int length)
        {
            Index = index;
            Length = length;
        }

        public int Index { get; private set; }
        public int Length { get; private set; }
    }
}