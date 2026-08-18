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
                // Only lazy-load from disk when there is a file to read. A factory-built
                // file (empty/non-existent FilePath) with Results set in memory returns
                // them directly instead of re-invoking LoadResults on a path that isn't
                // there -- which would crash, and is why callers no longer need a second
                // in-memory record store. File.Exists("" or null) is false (no throw).
                if (!_results.Any() && File.Exists(FilePath))
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

        #region Operators

        public TResult this[int index]
        {
            get => Results[index];
            set => Results[index] = value;
        }

        public static ResultFile<TResult> operator +(ResultFile<TResult> thisFile, TResult resultToAdd)
        {
            thisFile.Results.Add(resultToAdd);
            return thisFile;
        }

        public static ResultFile<TResult> operator +(ResultFile<TResult> thisFile, IEnumerable<TResult> resultsToAdd)
        {
            thisFile.Results.AddRange(resultsToAdd);
            return thisFile;
        }

        public static ResultFile<TResult> operator +(ResultFile<TResult> thisFile, ResultFile<TResult> fileToAdd)
        {
            thisFile.Results.AddRange(fileToAdd.Results);
            return thisFile;
        }

        public static ResultFile<TResult> operator -(ResultFile<TResult> thisFile, TResult resultToRemove)
        {
            thisFile.Results.Remove(resultToRemove);
            return thisFile;
        }

        public static ResultFile<TResult> operator -(ResultFile<TResult> thisFile, IEnumerable<TResult> resultsToRemove)
        {
            foreach (var result in resultsToRemove)
                thisFile.Results.Remove(result);
            return thisFile;
        }

        public static ResultFile<TResult> operator -(ResultFile<TResult> thisFile, ResultFile<TResult> fileToRemove)
        {
            foreach (var result in fileToRemove.Results)
                thisFile.Results.Remove(result);
            return thisFile;
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
