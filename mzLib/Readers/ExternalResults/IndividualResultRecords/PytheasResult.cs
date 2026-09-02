using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MzLibUtil;

namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// A single match line of a Pytheas match-output file. The file is not a delimited table
/// CsvHelper would read: each match line is space-delimited with a fixed 16-column prefix
/// followed by a variable number of MS2 match tokens and one final SCORE token, and the lines
/// are grouped under PRECURSOR_ION headers. Reading therefore has to split each line by hand.
/// </summary>
public class PytheasResult
{
    /// <summary>The m/z printed on the PRECURSOR_ION line this match belongs to.</summary>
    public string PrecursorIon { get; set; }

    /// <summary>The original match line, kept verbatim so that writing round-trips the file.</summary>
    public string RawLine { get; set; }

    /// <summary>m/z(meas): measured precursor m/z, minus any trailing '*' flag.</summary>
    public double MeasuredMz { get; set; }

    /// <summary>RT: retention time in minutes.</summary>
    public double RetentionTime { get; set; }

    /// <summary>m/z(theo): theoretical m/z of the matched oligomer.</summary>
    public double TheoreticalMz { get; set; }

    /// <summary>offset: mass error in parts per million.</summary>
    public double OffsetPpm { get; set; }

    /// <summary>Sp: Pytheas' match score.</summary>
    public double SpScore { get; set; }

    /// <summary>dSp: difference between the top and this match's score.</summary>
    public double DspScore { get; set; }

    /// <summary>rank: position of this match among the ones for its precursor ion.</summary>
    public int Rank { get; set; }

    /// <summary>#MS2_matches: number of MS2 matches reported for this precursor.</summary>
    public int Ms2MatchCount { get; set; }

    /// <summary>isotope: light or heavy.</summary>
    public string Isotope { get; set; }

    /// <summary>length: number of nucleotides in the matched oligomer.</summary>
    public int Length { get; set; }

    /// <summary>charge: precursor charge (negative).</summary>
    public int Charge { get; set; }

    /// <summary>sequence: oligomer sequence, e.g. "CCG".</summary>
    public string Sequence { get; set; }

    /// <summary>sequence_mod: modifications or ambiguity, e.g. "-", "CC[mG]CG" or "UUG|UCG".</summary>
    public string SequenceModification { get; set; }

    /// <summary>5'-end: terminal group at the 5' end, e.g. "OH".</summary>
    public string FivePrimeEnd { get; set; }

    /// <summary>3'-end: terminal group at the 3' end, e.g. "P".</summary>
    public string ThreePrimeEnd { get; set; }

    /// <summary>molecule_location: semicolon-delimited genomic locations, or "decoy".</summary>
    public string MoleculeLocation { get; set; }

    /// <summary>
    /// The raw MS2 match tokens, joined with single spaces in the order they appear, e.g.
    /// "892.167114(0.8ppm)[680]:892.166417[M-P(-1)] 362.051178(2.8ppm)[100]:362.050174[y1(-1)]".
    /// Kept as text rather than parsed so the line round-trips unchanged.
    /// </summary>
    public string Ms2Matches { get; set; }

    /// <summary>
    /// The MS2 <see cref="Ms2Matches"/> tokens parsed into typed <see cref="PytheasMatchedIon"/> records.
    /// Populated only when the tokens are requested; reading stays a single pass over the raw line.
    /// </summary>
    public List<PytheasMatchedIon> MatchedIons => _matchedIons ??= ParseMatchedIons();

    private List<PytheasMatchedIon>? _matchedIons;

    private List<PytheasMatchedIon> ParseMatchedIons()
    {
        if (string.IsNullOrWhiteSpace(Ms2Matches))
            return new List<PytheasMatchedIon>();
        return Ms2Matches.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(PytheasMatchedIon.Parse)
            .ToList();
    }

    /// <summary>The final SCORE token verbatim, e.g. "SCORE=0.239(sumI=227;n=5;...)".</summary>
    public string DetailedScore { get; set; }

    /// <summary>Numeric value extracted from the leading "SCORE=" part of <see cref="DetailedScore"/>.</summary>
    public double Score { get; set; }

    #region Interpreted Fields

    /// <summary>A match is a decoy when it has no genomic location, which Pytheas writes as "decoy".</summary>
    public bool IsDecoy => string.Equals(MoleculeLocation, "decoy", StringComparison.OrdinalIgnoreCase);

    /// <summary>The unmodified oligomer sequence, as printed in the sequence column.</summary>
    public string BaseSequence => Sequence;

    /// <summary>
    /// The sequence with modification or ambiguity annotations. Pytheas puts the plain sequence in its
    /// sequence column and the annotated form in sequence_mod; unmodified matches have a "-" placeholder.
    /// </summary>
    public string FullSequence => string.IsNullOrWhiteSpace(SequenceModification) || SequenceModification == "-"
        ? Sequence
        : SequenceModification;

    /// <summary>The genomic location (or "decoy") of the matched oligomer.</summary>
    public string Accession => MoleculeLocation;

    #endregion

    /// <summary>
    /// Parses one match line. The first 16 space-delimited tokens are the fixed columns
    /// documented by Pytheas' #MATCHES_HEADER line; the tokens up to the second-to-last are the
    /// MS2 matches; the last token is the SCORE.
    /// </summary>
    public static PytheasResult Parse(string line, string precursorIon)
    {
        string[] tokens = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
        if (tokens.Length < 17)
            throw new MzLibException($"Pytheas match line has too few columns: {line}");

        var result = new PytheasResult
        {
            PrecursorIon = precursorIon,
            RawLine = line,
            MeasuredMz = ParseDouble(tokens[0].TrimEnd('*')),
            RetentionTime = ParseDouble(ValueAfter(tokens[1], "RT=")),
            TheoreticalMz = ParseDouble(ValueAfter(tokens[2], "TH_MATCH=")),
            OffsetPpm = ParseDouble(TrimPpm(tokens[3])),
            SpScore = ParseDouble(ValueAfter(tokens[4], "Sp=")),
            DspScore = ParseDouble(ValueAfter(tokens[5], "dSp=")),
            Rank = ParseInt(ValueAfter(tokens[6], "rank=")),
            Ms2MatchCount = ParseInt(ValueAfter(tokens[7], "#MS2=")),
            Isotope = tokens[8],
            Length = ParseInt(tokens[9]),
            Charge = ParseInt(tokens[10]),
            Sequence = tokens[11],
            SequenceModification = tokens[12],
            FivePrimeEnd = tokens[13],
            ThreePrimeEnd = tokens[14],
            MoleculeLocation = tokens[15],
        };

        result.Ms2Matches = string.Join(" ", tokens, 16, tokens.Length - 17);
        result.DetailedScore = tokens[tokens.Length - 1];
        result.Score = ParseScore(result.DetailedScore);
        return result;
    }

    private static string TrimPpm(string token) =>
        token.EndsWith("ppm", StringComparison.OrdinalIgnoreCase) ? token.Substring(0, token.Length - 3) : token;

    private static string ValueAfter(string token, string prefix) =>
        token.StartsWith(prefix, StringComparison.Ordinal) ? token.Substring(prefix.Length) : token;

    private static double ParseScore(string detailedScore)
    {
        string value = ValueAfter(detailedScore, "SCORE=");
        if (value.Length == 0) return 0;
        int openParen = value.IndexOf('(');
        if (openParen >= 0) value = value.Substring(0, openParen);
        return ParseDouble(value);
    }

    private static double ParseDouble(string value) =>
        double.Parse(value, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static int ParseInt(string value) =>
        int.Parse(value, CultureInfo.InvariantCulture);
}