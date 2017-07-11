using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics
{
    public class DisulfideBond
    {
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

        #region Public Constructor

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

        #endregion Public Constructor
    }
}
