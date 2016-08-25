namespace Proteomics
{
    public class DigestionPoint
    {
        public int Index { get; private set; }
        public int Length { get; private set; }

        public DigestionPoint(int index, int length)
        {
            Index = index;
            Length = length;
        }
    }
}