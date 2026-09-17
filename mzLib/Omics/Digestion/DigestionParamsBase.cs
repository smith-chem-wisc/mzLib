using Omics.Fragmentation;

namespace Omics.Digestion
{
    /// <summary>
    /// Shared effective-agent state for proteolytic and nucleolytic digestion parameters.
    /// </summary>
    public abstract class DigestionParamsBase
    {
        public string DigestionAgentName => SpecificDigestionAgent.Name;
        public DigestionAgent DigestionAgent { get; protected set; }
        public DigestionAgent SpecificDigestionAgent { get; protected set; }
        public int MaxMissedCleavages { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MaxMods { get; set; }

        public FragmentationTerminus FragmentationTerminus { get; protected set; }
        public CleavageSpecificity SearchModeType { get; protected set; } = CleavageSpecificity.Full;

        public DigestionParamsBase(DigestionAgent digestionAgent, int maxMissedCleavages = 0, int minLength = 3, int maxLength = int.MaxValue, int maxModificationIsoforms = 1024, int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full)
        {
            MaxMissedCleavages = maxMissedCleavages;
            MinLength = minLength;
            MaxLength = maxLength;
            MaxMods = maxMods;
            MaxModificationIsoforms = maxModificationIsoforms;
            FragmentationTerminus = fragmentationTerminus;
            SearchModeType = searchModeType;

            DigestionAgent = digestionAgent;
            if (SearchModeType == CleavageSpecificity.None)
            {
                SpecificDigestionAgent = digestionAgent;
                DigestionAgent = GetSingleTerminusAgent(FragmentationTerminus);
            }
            else
            {
                SpecificDigestionAgent = digestionAgent;
            }
        }

        protected abstract DigestionAgent GetSingleTerminusAgent(FragmentationTerminus terminus);
    }
}
