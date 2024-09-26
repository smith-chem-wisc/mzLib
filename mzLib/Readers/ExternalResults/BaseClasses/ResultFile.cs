using System;
using System.Collections;

namespace Readers
{
    /// <summary>
    /// Generic abstract class for all types of result files in MzLib
    /// </summary>
    /// <typeparam name="TResult">type of result being contained</typeparam>
    public abstract class ResultFile<TResult> : IResultFile, IEquatable<ResultFile<TResult>>, IEnumerable<TResult> 
    {
        #region Base Properties
        public string FilePath { get; set; }

        private List<TResult> _results;
        public List<TResult> Results
        {
            get
            {
                if (!_results.Any())
                    LoadResults();
                return _results;
            }
            set => _results = value;
        }

        #endregion

        #region Constructors

        protected ResultFile(string filePath, Software software = Software.Unspecified)
        {
            FilePath = filePath;
            _results = new List<TResult>();
            Software = software;
        }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        protected internal ResultFile()
        {
            FilePath = "";
            _results = new();
            Software = Software.Unspecified;
        }

        #endregion

        #region Abstract Members

        public abstract SupportedFileType FileType { get; }
        public abstract Software Software { get; set; } 

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public abstract void LoadResults();

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public abstract void WriteResults(string outputPath);

        #endregion

        #region Methods

        /// <summary>
        /// Determines whether the specified file type and extension of filepath align
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public bool CanRead(string filePath)
        {
            return filePath.EndsWith(FileType.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase);
        }

        #endregion

        #region Interface Implementations

        public IEnumerator<TResult> GetEnumerator() => Results.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        public bool Equals(ResultFile<TResult>? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return GetHashCode() == other.GetHashCode();
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((ResultFile<TResult>)obj);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(FilePath);
        }

        #endregion



    }
}
