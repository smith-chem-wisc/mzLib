using System.Collections.Generic;

namespace UsefulProteomicsDatabases.Transcriptomics;

/// <summary>
/// Intermediate Data Transfer Object (DTO) for Modomics entries.
/// </summary>
public class ModomicsDto
{
    public int Id { get; set; }
    public string Abbrev { get; set; }
    public string Formula { get; set; }
    public string LcElutionComment { get; set; }
    public string LcElutionTime { get; set; }
    public double MassAvg { get; set; }
    public double MassMonoiso { get; set; }
    public double? MassProt { get; set; }
    public string Name { get; set; }
    public string ProductIons { get; set; }
    public List<string> ReferenceMoiety { get; set; }
    public string ShortName { get; set; }
    public string Smile { get; set; }
    public string MoietyType { get; set; }
}