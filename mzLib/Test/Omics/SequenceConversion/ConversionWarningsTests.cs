using NUnit.Framework;
using Omics.SequenceConversion;
using System.Linq;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ConversionWarningsTests
{
    [Test]
    public void NewInstance_IsCleanAndHasNoIssues()
    {
        var warnings = new ConversionWarnings();

        Assert.That(warnings.IsClean, Is.True);
        Assert.That(warnings.HasWarnings, Is.False);
        Assert.That(warnings.HasErrors, Is.False);
        Assert.That(warnings.HasIncompatibleItems, Is.False);
        Assert.That(warnings.HasFatalError, Is.False);
        Assert.That(warnings.ToString(), Is.EqualTo("No issues"));
    }

    [Test]
    public void AddMethods_IgnoreWhitespaceAndTrackCounts()
    {
        var warnings = new ConversionWarnings();

        warnings.AddWarning("warn");
        warnings.AddWarning(" ");
        warnings.AddError("error");
        warnings.AddError("");
        warnings.AddIncompatibleItem("item");
        warnings.AddIncompatibleItem("\t");

        Assert.That(warnings.Warnings.Count, Is.EqualTo(1));
        Assert.That(warnings.Errors.Count, Is.EqualTo(1));
        Assert.That(warnings.IncompatibleItems.Count, Is.EqualTo(1));
    }

    [Test]
    public void Merge_PreservesExistingFailureReasonAndAppendsIssues()
    {
        var left = new ConversionWarnings();
        left.AddWarning("left-warning");
        left.SetFailure(ConversionFailureReason.IncompatibleModifications, "left-error");

        var right = new ConversionWarnings();
        right.AddWarning("right-warning");
        right.AddError("right-error");
        right.AddIncompatibleItem("right-item");
        right.SetFailure(ConversionFailureReason.InvalidSequence, "right-fatal");

        left.Merge(right);

        Assert.That(left.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
        Assert.That(left.Warnings.Count, Is.EqualTo(2));
        Assert.That(left.Errors.Count, Is.EqualTo(3));
        Assert.That(left.IncompatibleItems.Single(), Is.EqualTo("right-item"));
    }

    [Test]
    public void ToException_UsesFailureReasonAndCapturedLists()
    {
        var warnings = new ConversionWarnings();
        warnings.AddWarning("warn");
        warnings.AddIncompatibleItem("item");
        warnings.SetFailure(ConversionFailureReason.IncompatibleModifications, "error");

        var ex = warnings.ToException("boom");

        Assert.That(ex.Message, Does.Contain("boom"));
        Assert.That(ex.FailureReason, Is.EqualTo(ConversionFailureReason.IncompatibleModifications));
        Assert.That(ex.IncompatibleItems, Is.Not.Null);
        Assert.That(ex.IncompatibleItems!.Single(), Is.EqualTo("item"));
        Assert.That(ex.Warnings, Is.Not.Null);
        Assert.That(ex.Warnings!.Single(), Is.EqualTo("warn"));
    }

    [Test]
    public void Clear_RemovesAllState()
    {
        var warnings = new ConversionWarnings();
        warnings.AddWarning("warn");
        warnings.AddError("error");
        warnings.AddIncompatibleItem("item");
        warnings.SetFailure(ConversionFailureReason.InvalidSequence, "fatal");

        warnings.Clear();

        Assert.That(warnings.IsClean, Is.True);
        Assert.That(warnings.Warnings, Is.Empty);
        Assert.That(warnings.Errors, Is.Empty);
        Assert.That(warnings.IncompatibleItems, Is.Empty);
        Assert.That(warnings.FailureReason, Is.Null);
    }
}
