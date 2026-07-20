using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Omics.BioPolymerGroup;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Renders groups through the production schema and writer.
    ///
    /// Groups no longer format themselves, so tests go through the same path real output does.
    /// The single-group overloads build a schema from just that group, which matches how these
    /// tests were originally written; <see cref="RowInDataset"/> exists for the cases that need a
    /// row rendered against a schema derived from several groups.
    /// </summary>
    [ExcludeFromCodeCoverage]
    internal static class GroupTsv
    {
        /// <summary>Header line for a file containing the given groups.</summary>
        public static string Header(params BioPolymerGroup[] groups) =>
            TsvWriter.HeaderLine(BioPolymerGroupTsvSchema.For(groups));

        /// <summary>Row for a group written on its own.</summary>
        public static string Row(BioPolymerGroup group) =>
            TsvWriter.RowLine(BioPolymerGroupTsvSchema.For([group]), group);

        /// <summary>Row for one group rendered against the schema of the whole dataset.</summary>
        public static string RowInDataset(BioPolymerGroup group, IReadOnlyCollection<BioPolymerGroup> dataset) =>
            TsvWriter.RowLine(BioPolymerGroupTsvSchema.For(dataset), group);

        /// <summary>Header line for a file containing the given digestion-product groups.</summary>
        public static string PeptideHeader(params BioPolymerWithSetModsGroup[] groups) =>
            TsvWriter.HeaderLine(BioPolymerWithSetModsGroupTsvSchema.For(groups));

        /// <summary>Row for a digestion-product group written on its own.</summary>
        public static string PeptideRow(BioPolymerWithSetModsGroup group) =>
            TsvWriter.RowLine(BioPolymerWithSetModsGroupTsvSchema.For([group]), group);

        /// <summary>Row for one digestion-product group rendered against the whole dataset's schema.</summary>
        public static string PeptideRowInDataset(
            BioPolymerWithSetModsGroup group, IReadOnlyCollection<BioPolymerWithSetModsGroup> dataset) =>
            TsvWriter.RowLine(BioPolymerWithSetModsGroupTsvSchema.For(dataset), group);
    }
}
