#nullable enable
using System;
using System.Collections.Generic;
using System.Threading;
using MassSpectrometry;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace PredictionClients.SpectralLibraryGeneration
{
    public sealed record MslLibraryBuildParameters(
        DigestionParams DigestionParams,
        List<Modification> FixedModifications,
        List<Modification> VariableModifications,
        IReadOnlyList<int> PrecursorCharges,
        int Nce,
        DissociationType DissociationType,
        DecoyType DecoyType = DecoyType.Reverse,
        int PredictionChunkSize = 50_000);

    public sealed record MslLibraryBuildResult(int TargetPrecursors, int DecoyPrecursors, int NotPredicted, int NoFragments,
        int DecoysDroppedAsTargetSequences);

    public sealed class MslLibraryBuilder
    {
        public MslLibraryBuilder(FragmentIntensityModel intensityModel, RetentionTimeModel irtModel)
        {
        }

        public Func<IReadOnlyList<Protein>, IReadOnlyList<Protein>>? EntrapmentProteins { get; init; }

        public List<MslLibraryEntry> Build(IReadOnlyList<Protein> targets, MslLibraryBuildParameters parameters,
            out MslLibraryBuildResult result, CancellationToken cancellationToken = default)
        {
            throw new NotImplementedException();
        }
    }
}
