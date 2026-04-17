using Omics.Modifications;

namespace Omics.SequenceConversion;

public class KnownModsLookup : ModificationLookupBase
{
    public KnownModsLookup(Dictionary<string, Modification> knownMods, double massTolerance = 0.001)
        : base(knownMods?.Values, massTolerance)
    {
    }

    public override string Name => "Known Mods";
}
