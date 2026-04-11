namespace Readers;

public class ProSightPdProteinRecord
{
    public int Id { get; set; }
    public string? Accession { get; set; }
    public int IsoformCount { get; set; }
    public int IsoformsWithCharacterizedProteoformsCount { get; set; }
    public int ProteoformCount { get; set; }
    public int ProteoformCharacterizationConfidence { get; set; }
    public int CharacterizedProteoformsCount { get; set; }
    public int PrSMCount { get; set; }
    public int ExcludedBy { get; set; }
    public string? Description { get; set; }
    public double Qvalue { get; set; }
}
