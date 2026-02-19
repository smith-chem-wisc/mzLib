using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;

namespace Test.Omics;

public class MockBioPolymer : IBioPolymer
{
    public string BaseSequence { get; }
    public string Accession { get; }
    public string Organism { get; }
    public string Name { get; }
    public string FullName { get; }
    public List<Tuple<string, string>> GeneNames { get; }
    public bool IsDecoy { get; }
    public bool IsContaminant { get; }
    public string DatabaseFilePath { get; } = "";
    public int Length => BaseSequence.Length;
    public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; } = new Dictionary<int, List<Modification>>();
    public string SampleNameForVariants { get; set; } = "";
    public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; } = new Dictionary<int, List<Modification>>();
    public IBioPolymer ConsensusVariant => this;
    public List<SequenceVariation> AppliedSequenceVariations { get; } = new();
    public List<SequenceVariation> SequenceVariations { get; } = new();
    public List<TruncationProduct> TruncationProducts { get; } = new();

    public MockBioPolymer(string sequence, string accession) : this(sequence, accession, false, false){}

    public MockBioPolymer(string sequence, string accession, bool isDecoy = false, bool isContaminant = false) : this(sequence, accession, "", "", "", null, isDecoy, isContaminant){}

    public MockBioPolymer(string sequence, string accession,
        string organism = "", string name = "", string fullName = "",
        List<Tuple<string, string>> geneNames = null,
        bool isDecoy = false, bool isContaminant = false)
    {
        BaseSequence = sequence;
        Accession = accession;
        Organism = organism;
        Name = name;
        FullName = fullName;
        GeneNames = geneNames ?? new List<Tuple<string, string>>();
        IsDecoy = isDecoy;
        IsContaminant = isContaminant;
    }

    public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
        List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
        List<SilacLabel>? silacLabels = null, (SilacLabel, SilacLabel)? turnoverLabels = null,
        bool topDownTruncationSearch = false) => Enumerable.Empty<IBioPolymerWithSetMods>();

    public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods)
        => new MockBioPolymer(newBaseSequence, Accession, Organism, Name, FullName, GeneNames, IsDecoy, IsContaminant);

    public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
        IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts,
        IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
        where TBioPolymerType : IHasSequenceVariants => original;

    public bool Equals(IBioPolymer? other) => other != null && Accession == other.Accession && BaseSequence == other.BaseSequence;
    public override bool Equals(object? obj) => obj is IBioPolymer other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
}