using System;
using System.Collections.Generic;
using System.IO;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>One pathway of a Reactome release.</summary>
    public sealed record ReactomePathwayRow(string PathwayId, string Name, string SpeciesName, string HumanEquivalentId);

    /// <summary>One parent-child pair of a Reactome release's hierarchy.</summary>
    public sealed record ReactomeHierarchyRow(string ParentId, string ChildId);

    /// <summary>One pair of a Reactome release's ancestor closure.</summary>
    public sealed record ReactomeAncestorRow(string PathwayId, string AncestorPathwayId);

    // RED: surface only, so the failing tests compile.
    public static class ReactomeReleaseTsv
    {
        public const int FormatVersion = 1;
        public static void WritePathways(TextWriter output, ReactomeRelease release) => throw new NotImplementedException();
        public static void WriteHierarchy(TextWriter output, ReactomeRelease release) => throw new NotImplementedException();
        public static void WriteAncestors(TextWriter output, ReactomeRelease release) => throw new NotImplementedException();
        public static AnnotationTable<ReactomePathwayRow> ReadPathways(string path) => throw new NotImplementedException();
        public static AnnotationTable<ReactomeHierarchyRow> ReadHierarchy(string path) => throw new NotImplementedException();
        public static AnnotationTable<ReactomeAncestorRow> ReadAncestors(string path) => throw new NotImplementedException();
    }
}
