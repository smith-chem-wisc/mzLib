using Omics.Fragmentation;

namespace Omics.Digestion
{
    public interface IDigestionParams 
    {
        int MaxMissedCleavages { get; set; }
        int MinLength { get; set; }
        int MaxLength { get; set; }
        int MaxModificationIsoforms { get; set; }
        int MaxMods { get; set; }
        DigestionAgent DigestionAgent { get; }
        FragmentationTerminus FragmentationTerminus { get; }
        CleavageSpecificity SearchModeType { get; }

        IDigestionParams Clone();
    }
}


//(CommonParameters.DigestionParams as DigestionParams)!
//if (CommonParameters.DigestionParams is DigestionParams digestion)
//    _ = ProseCreatedWhileRunning.Append("initiator methionine behavior = " + digestion.InitiatorMethionineBehavior + "; ");