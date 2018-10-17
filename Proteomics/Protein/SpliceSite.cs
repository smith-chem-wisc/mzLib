namespace Proteomics
{
    public class SpliceSite
    {
        public SpliceSite(int oneBasedBegin, int oneBasedEnd, string description)
        {
            OneBasedBeginPosition = oneBasedBegin;
            OneBasedEndPosition = oneBasedEnd;
            Description = description;
        }

        public SpliceSite(int oneBasedPosition, string description)
            : this(oneBasedPosition, oneBasedPosition, description)
        {
        }

        public int OneBasedBeginPosition { get; }
        public int OneBasedEndPosition { get; }
        public string Description { get; }
    }
}