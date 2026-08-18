namespace Omics.SequenceConversion;

/// <summary>
/// Accumulates warnings and errors during sequence conversion operations.
/// Passed as a nullable parameter to parsing/serialization methods and built up through sub-steps.
/// </summary>
public class ConversionWarnings
{
    private readonly List<string> _warnings = new();
    private readonly List<string> _errors = new();
    private readonly List<string> _incompatibleItems = new();

    /// <summary>
    /// Gets the list of warning messages accumulated during conversion.
    /// Warnings indicate non-fatal issues that didn't prevent conversion.
    /// </summary>
    public IReadOnlyList<string> Warnings => _warnings;

    /// <summary>
    /// Gets the list of error messages accumulated during conversion.
    /// Errors indicate issues that may have affected the conversion result.
    /// </summary>
    public IReadOnlyList<string> Errors => _errors;

    /// <summary>
    /// Gets the list of incompatible items encountered during conversion.
    /// These are specific modifications or elements that couldn't be converted.
    /// </summary>
    public IReadOnlyList<string> IncompatibleItems => _incompatibleItems;

    /// <summary>
    /// Gets the failure reason if a fatal error occurred.
    /// Null if no fatal error has been recorded.
    /// </summary>
    public ConversionFailureReason? FailureReason { get; private set; }

    /// <summary>
    /// Returns true if any warnings have been recorded.
    /// </summary>
    public bool HasWarnings => _warnings.Count > 0;

    /// <summary>
    /// Returns true if any errors have been recorded.
    /// </summary>
    public bool HasErrors => _errors.Count > 0;

    /// <summary>
    /// Returns true if any incompatible items have been recorded.
    /// </summary>
    public bool HasIncompatibleItems => _incompatibleItems.Count > 0;

    /// <summary>
    /// Returns true if a fatal failure reason has been set.
    /// </summary>
    public bool HasFatalError => FailureReason.HasValue;

    /// <summary>
    /// Returns true if the conversion completed without any issues.
    /// </summary>
    public bool IsClean => !HasWarnings && !HasErrors && !HasIncompatibleItems && !HasFatalError;

    /// <summary>
    /// Adds a warning message.
    /// </summary>
    public void AddWarning(string message)
    {
        if (!string.IsNullOrWhiteSpace(message))
            _warnings.Add(message);
    }

    /// <summary>
    /// Adds an error message.
    /// </summary>
    public void AddError(string message)
    {
        if (!string.IsNullOrWhiteSpace(message))
            _errors.Add(message);
    }

    /// <summary>
    /// Adds an incompatible item description.
    /// </summary>
    public void AddIncompatibleItem(string item)
    {
        if (!string.IsNullOrWhiteSpace(item))
            _incompatibleItems.Add(item);
    }

    /// <summary>
    /// Sets a fatal failure reason. This indicates the conversion could not complete.
    /// </summary>
    public void SetFailure(ConversionFailureReason reason, string? errorMessage = null)
    {
        FailureReason = reason;
        if (!string.IsNullOrWhiteSpace(errorMessage))
            _errors.Add(errorMessage);
    }

    /// <summary>
    /// Merges another ConversionWarnings instance into this one.
    /// Useful for aggregating warnings from sub-operations.
    /// </summary>
    public void Merge(ConversionWarnings? other)
    {
        if (other == null) return;

        _warnings.AddRange(other._warnings);
        _errors.AddRange(other._errors);
        _incompatibleItems.AddRange(other._incompatibleItems);

        // Only overwrite failure reason if we don't have one and the other does
        if (!FailureReason.HasValue && other.FailureReason.HasValue)
            FailureReason = other.FailureReason;
    }

    /// <summary>
    /// Clears all accumulated warnings, errors, and failure state.
    /// </summary>
    public void Clear()
    {
        _warnings.Clear();
        _errors.Clear();
        _incompatibleItems.Clear();
        FailureReason = null;
    }

    /// <summary>
    /// Creates a SequenceConversionException from the current state.
    /// </summary>
    public SequenceConversionException ToException(string message)
    {
        return new SequenceConversionException(
            message,
            FailureReason ?? ConversionFailureReason.UnknownFormat,
            _incompatibleItems.Count > 0 ? _incompatibleItems : null,
            _warnings.Count > 0 ? _warnings : null);
    }

    /// <summary>
    /// Returns a summary string of all issues.
    /// </summary>
    public override string ToString()
    {
        var parts = new List<string>();

        if (HasFatalError)
            parts.Add($"Fatal: {FailureReason}");
        if (HasErrors)
            parts.Add($"Errors: {_errors.Count}");
        if (HasWarnings)
            parts.Add($"Warnings: {_warnings.Count}");
        if (HasIncompatibleItems)
            parts.Add($"Incompatible: {_incompatibleItems.Count}");

        return parts.Count > 0 ? string.Join(", ", parts) : "No issues";
    }
}
