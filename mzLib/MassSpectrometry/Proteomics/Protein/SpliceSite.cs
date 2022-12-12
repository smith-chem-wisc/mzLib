namespace Proteomics
{
    public class SpliceSite
    {
        public SpliceSite(int oneBasedBegin, int oneBasedEnd, string description)
        {
            OneBasedBeginPosition = oneBasedBegin;
            OneBasedEndPosition = oneBasedEnd;
            Description = description ?? "";
        }

        public SpliceSite(int oneBasedPosition, string description)
            : this(oneBasedPosition, oneBasedPosition, description)
        {
        }

        public int OneBasedBeginPosition { get; }
        public int OneBasedEndPosition { get; }
        public string Description { get; }

        public override bool Equals(object obj)
        {
            SpliceSite s = obj as SpliceSite;
            return s != null
                && s.OneBasedBeginPosition == OneBasedBeginPosition
                && s.OneBasedEndPosition == OneBasedEndPosition
                && s.Description == Description;
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode() 
                ^ OneBasedEndPosition.GetHashCode() 
                ^ Description.GetHashCode(); // null handled in constructor
        }
    }
}