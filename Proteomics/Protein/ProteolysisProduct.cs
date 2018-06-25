namespace Proteomics
{
    public class ProteolysisProduct
    {
        public ProteolysisProduct(int? oneBasedBeginPosition, int? oneBasedEndPosition, string type)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            Type = type ?? "";
        }

        public int? OneBasedBeginPosition { get; }
        public int? OneBasedEndPosition { get; }
        public string Type { get; }

        public override bool Equals(object obj)
        {
            return obj as ProteolysisProduct != null &&
                (obj as ProteolysisProduct).OneBasedBeginPosition == OneBasedBeginPosition &&
                (obj as ProteolysisProduct).OneBasedEndPosition == OneBasedEndPosition &&
                (obj as ProteolysisProduct).Type == Type;
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode() ^ OneBasedEndPosition.GetHashCode() ^ Type.GetHashCode();
        }
    }
}