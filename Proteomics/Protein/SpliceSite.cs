namespace Proteomics
{
    public class SpliceSite
    {
        public SpliceSite(int oneBasedBegin, int oneBasedEnd, string description)
        {
            OneBasedBeginPosition = oneBasedBegin;
            OneBasedEndPosition = oneBasedEnd;
            Description = new SpliceSiteDescription(description);
        }

        public SpliceSite(int oneBasedPosition, string description)
            : this(oneBasedPosition, oneBasedPosition, description)
        {
        }

        public int OneBasedBeginPosition { get; }
        public int OneBasedEndPosition { get; }
        public SpliceSiteDescription Description { get; }

        public override bool Equals(object obj)
        {
            SpliceSite s = obj as SpliceSite;
            return s != null
                && s.OneBasedBeginPosition == OneBasedBeginPosition
                && s.OneBasedEndPosition == OneBasedEndPosition
                && s.Description.Description == Description.Description;
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode() 
                ^ OneBasedEndPosition.GetHashCode() 
                ^ Description.GetHashCode(); // null handled in constructor
        }
    }
}