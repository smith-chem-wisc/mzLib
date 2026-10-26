using System;
using System.Collections.Generic;

namespace UsefulProteomicsDatabases.Pathways
{
    // RED: surface only, so the failing tests compile.
    public sealed class ReactomeRelease
    {
        public string Release { get; }
        public string PathwaysSha256 { get; }
        public string RelationSha256 { get; }
        public IReadOnlyList<ReactomePathwayRow> Pathways { get; } = Array.Empty<ReactomePathwayRow>();
        public IReadOnlyList<ReactomeHierarchyRow> Hierarchy { get; } = Array.Empty<ReactomeHierarchyRow>();
        public IReadOnlyList<string> Roots { get; } = Array.Empty<string>();
        public bool Contains(string pathwayId) => throw new NotImplementedException();
        public IReadOnlySet<string> Ancestors(string pathwayId) => throw new NotImplementedException();
        public string HumanEquivalentId(string pathwayId) => throw new NotImplementedException();
        public static ReactomeRelease Load(string pathwaysPath, string relationPath, string release) => throw new NotImplementedException();
    }
}
