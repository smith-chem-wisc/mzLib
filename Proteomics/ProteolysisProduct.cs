namespace Proteomics
{
    public class ProteolysisProduct
    {

        #region Public Constructors

        public ProteolysisProduct(int? oneBasedBeginPosition, int? oneBasedEndPosition, string type)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            Type = type;
        }

        #endregion Public Constructors

        #region Public Properties

        public int? OneBasedBeginPosition { get; }
        public int? OneBasedEndPosition { get; }
        public string Type { get; }

        #endregion Public Properties

        #region Public Methods

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

        #endregion Public Methods

    }
}