﻿namespace Proteomics
{
    public class DisulfideBond
    {
        #region Public Constructors

        public DisulfideBond(int OneBasedBeginPosition, int OneBasedEndPosition, string Description)
        {
            this.OneBasedBeginPosition = OneBasedBeginPosition;
            this.OneBasedEndPosition = OneBasedEndPosition;
            this.Description = Description;
        }

        /// For interchain disulfide bonds, sets begin and end to the same position.
        public DisulfideBond(int OneBasedPosition, string Description)
            : this(OneBasedPosition, OneBasedPosition, Description)
        { }

        #endregion Public Constructors

        #region Public Properties

        /// <summary>
        /// Beginning position of disulfide bond
        /// </summary>
        public int OneBasedBeginPosition { get; set; }

        /// <summary>
        /// End position of disulfide bond
        /// </summary>
        public int OneBasedEndPosition { get; set; }

        /// <summary>
        /// Description of this variation (optional)
        /// </summary>
        public string Description { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override bool Equals(object obj)
        {
            DisulfideBond bond = obj as DisulfideBond;
            return bond != null
                && bond.OneBasedBeginPosition == OneBasedBeginPosition
                && bond.OneBasedEndPosition == OneBasedEndPosition
                && bond.Description == Description;
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition ^ OneBasedEndPosition ^ Description.GetHashCode();
        }

        #endregion Public Methods
    }
}