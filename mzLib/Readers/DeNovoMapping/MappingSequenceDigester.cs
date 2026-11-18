using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Readers.DeNovoMapping;

/// <summary>
/// Receives the peptides of a digestion as windows of a normalized sequence, so no peptide object or string is made per
/// peptide.
/// </summary>
internal interface IMappingWindowSink
{
    /// <param name="normalizedSource">The I/L-normalized protein, or one normalized peptide.</param>
    void OnWindow(string normalizedSource, int start, int length);
}

/// <summary>
/// Enumerates the unmodified peptides <see cref="IBioPolymer.Digest"/> produces, normalized as
/// <see cref="DeNovoPrediction.NormalizeSequence"/> does, as windows of the normalized protein. For plain fully specific
/// digestion it walks the cleavage sites and truncation products exactly as <c>DigestionAgent.FullDigestion</c> does,
/// without building a <see cref="PeptideWithSetModifications"/> per peptide. Anything else (semi or non-specific
/// digestion, glycopeptide filters, an active <see cref="CleavageBlockingPolicy"/>, a truncation product outside the
/// protein, a protein that is not a <see cref="Protein"/>) goes through <see cref="IBioPolymer.Digest"/> itself. A
/// peptide may be reported more than once.
/// </summary>
internal sealed class MappingSequenceDigester
{
    private readonly DigestionParams _digestionParams;
    private readonly bool _replaceIWithL;
    private readonly bool _isFullDigestion;
    private readonly bool _policyCanActivate;
    private readonly List<Modification> _fixedModifications = new();
    private readonly List<Modification> _variableModifications = new();

    public MappingSequenceDigester(DigestionParams digestionParams, bool replaceIWithL)
    {
        _digestionParams = digestionParams;
        _replaceIWithL = replaceIWithL;

        // The lists Protein.Digest builds from the empty ones the mapper passes, so the policy check below sees the same.
        IBioPolymer.AddDigestionAgentModification(digestionParams.Protease, ref _fixedModifications, ref _variableModifications);

        // Protease.GetUnmodifiedPeptides takes FullDigestion exactly when neither the search mode nor the protease asks
        // for anything else; SearchModeType None has already swapped in a single-terminus protease.
        _isFullDigestion = digestionParams.SearchModeType == CleavageSpecificity.Full
                           && digestionParams.Protease.CleavageSpecificity == CleavageSpecificity.Full
                           && !digestionParams.KeepNGlycopeptide
                           && !digestionParams.KeepOGlycopeptide;

        // CleavageBlockingPolicy.For is inert unless one of these is on, so only then is it asked per protein.
        _policyCanActivate = digestionParams.RespectCleavageBlockingModifications || digestionParams.RespectCleavagePromotingModifications;
    }

    /// <summary>
    /// True when <paramref name="protein"/> is enumerated from its cleavage sites rather than through <see cref="IBioPolymer.Digest"/>.
    /// </summary>
    public bool EnumeratesWindows(IBioPolymer protein)
    {
        if (!_isFullDigestion || protein is not Protein { BaseSequence.Length: > 0 } concrete)
            return false;

        foreach (var truncation in concrete.TruncationProducts)
        {
            if (!IsInside(truncation, concrete.Length))
                return false;
        }

        return !_policyCanActivate || !CleavageBlockingPolicy.For(_digestionParams, _variableModifications, concrete, _fixedModifications).IsActive;
    }

    /// <summary>
    /// True when a truncation product's known ends lie inside the protein, in order, so FullDigestion's walk over it stays
    /// in bounds. Anything else is left to <see cref="IBioPolymer.Digest"/>, to behave exactly as it does.
    /// </summary>
    private static bool IsInside(TruncationProduct truncation, int proteinLength)
    {
        int? begin = truncation.OneBasedBeginPosition;
        int? end = truncation.OneBasedEndPosition;
        bool beginInside = begin is null || (begin >= 1 && begin <= proteinLength);
        bool endInside = end is null || (end >= 1 && end <= proteinLength);
        bool ordered = begin is null || end is null || begin <= end;
        return beginInside && endInside && ordered;
    }

    public void Digest(IBioPolymer protein, IMappingWindowSink sink)
    {
        if (!EnumeratesWindows(protein))
        {
            foreach (var peptide in protein.Digest(_digestionParams, new List<Modification>(), new List<Modification>()))
            {
                string normalizedPeptide = DeNovoPrediction.NormalizeSequence(peptide.BaseSequence, _replaceIWithL);
                sink.OnWindow(normalizedPeptide, 0, normalizedPeptide.Length);
            }
            return;
        }

        string sequence = protein.BaseSequence;
        string normalized = DeNovoPrediction.NormalizeSequence(sequence, _replaceIWithL);
        List<int> sites = _digestionParams.Protease.GetDigestionSiteIndices(sequence);

        // As Protease.GetUnmodifiedPeptides and DigestionAgent.FullDigestion, with sites as zero-based starts.
        var initiatorMethionine = _digestionParams.InitiatorMethionineBehavior;
        bool retainFirstResidue = initiatorMethionine != InitiatorMethionineBehavior.Cleave || sequence[0] != 'M';
        bool cleaveFirstResidue = initiatorMethionine != InitiatorMethionineBehavior.Retain && sequence[0] == 'M';
        int minLength = _digestionParams.MinLength;
        int maxLength = _digestionParams.MaxLength;

        int maxMissedCleavages = _digestionParams.MaxMissedCleavages;
        bool ValidLength(int length) => length >= minLength && length <= maxLength;

        // A truncation product without both ends adds no peptide in FullDigestion; EnumeratesWindows has checked the rest.
        var truncationProducts = protein.TruncationProducts.Count == 0
            ? []
            : protein.TruncationProducts
                .Where(t => t.OneBasedBeginPosition.HasValue && t.OneBasedEndPosition.HasValue)
                .Select(t => (Begin: t.OneBasedBeginPosition!.Value, End: t.OneBasedEndPosition!.Value))
                .ToList();

        for (int missedCleavages = 0; missedCleavages <= maxMissedCleavages; missedCleavages++)
        {
            for (int i = 0; i < sites.Count - missedCleavages - 1; i++)
            {
                int end = sites[i + missedCleavages + 1];
                if ((i != 0 || retainFirstResidue) && ValidLength(end - sites[i]))
                    sink.OnWindow(normalized, sites[i], end - sites[i]);

                if (i == 0 && cleaveFirstResidue && sites[1] != 1 && ValidLength(end - 1))
                    sink.OnWindow(normalized, 1, end - 1);
            }

            // Peptides that start or end at a truncation product's own terminus. One-based residues begin..end are the
            // window starting at begin - 1, of length end - begin + 1.
            foreach (var (begin, end) in truncationProducts)
            {
                if (begin == 1 && end == sequence.Length)
                    continue;

                int cleavageIndex = 0;
                while (sites[cleavageIndex] < begin)
                    cleavageIndex++;

                // From the truncation's start to the cleavage site missedCleavages sites on.
                if (cleavageIndex + missedCleavages < sites.Count)
                {
                    int lastResidue = sites[cleavageIndex + missedCleavages];
                    if (lastResidue <= end
                        && !sites.Contains(begin - 1)
                        && (begin != 1 || !cleaveFirstResidue)
                        && ValidLength(lastResidue - begin + 1))
                    {
                        sink.OnWindow(normalized, begin - 1, lastResidue - begin + 1);
                    }
                }

                while (sites[cleavageIndex] < end)
                    cleavageIndex++;

                // From the cleavage site missedCleavages sites back to the truncation's end.
                if (cleavageIndex - missedCleavages - 1 >= 0)
                {
                    int startSite = sites[cleavageIndex - missedCleavages - 1];
                    if (startSite + 1 >= begin && !sites.Contains(end) && ValidLength(end - startSite))
                        sink.OnWindow(normalized, startSite, end - startSite);
                }
            }
        }

        // Whole truncation products. FullDigestion checks end - begin, one less than the length, against the length limits.
        foreach (var (begin, end) in truncationProducts)
        {
            if ((begin == 1 && cleaveFirstResidue) || sites.Contains(begin - 1) || sites.Contains(end) || !ValidLength(end - begin))
                continue;

            int firstCleavage = 0;
            while (sites[firstCleavage] < begin)
                firstCleavage++;

            int lastCleavage = firstCleavage;
            while (sites[lastCleavage] < end)
                lastCleavage++;

            if (lastCleavage - firstCleavage < maxMissedCleavages)
                sink.OnWindow(normalized, begin - 1, end - begin + 1);
        }
    }
}
