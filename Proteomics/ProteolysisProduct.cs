using System;

namespace Proteomics
{
    [Serializable]
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

        public int? OneBasedBeginPosition { get; private set; }
        public int? OneBasedEndPosition { get; private set; }
        public string Type { get; private set; }

        #endregion Public Properties

    }
}