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

    }
}