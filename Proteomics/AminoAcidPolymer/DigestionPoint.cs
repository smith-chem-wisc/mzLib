namespace Proteomics.AminoAcidPolymer
{
    public class DigestionPointAndLength
    {
        #region Public Constructors

        public DigestionPointAndLength(int index, int length)
        {
            Index = index;
            Length = length;
        }

        #endregion Public Constructors

        #region Public Properties

        public int Index { get; private set; }
        public int Length { get; private set; }

        #endregion Public Properties
    }
}