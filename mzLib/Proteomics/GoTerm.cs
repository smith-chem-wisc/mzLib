namespace Proteomics
{
    public enum Aspect
    {
        molecularFunction,
        cellularComponent,
        biologicalProcess
    }

    public class GoTerm
    {

        #region Public Constructors

        public GoTerm(string id, string description, Aspect aspect)
        {
            Id = id;
            Description = description;
            Aspect = aspect;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Id { get; private set; }
        public string Description { get; private set; }
        public Aspect Aspect { get; private set; }

        #endregion Public Properties

    }
}