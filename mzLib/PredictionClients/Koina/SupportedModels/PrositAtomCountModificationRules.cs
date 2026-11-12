using System.Text.RegularExpressions;
using Chemistry;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace PredictionClients.Koina.SupportedModels;

/// <summary>
/// The modifications Koina's atom-count Prosit models (2024 PTMs_gl, 2025 22PTM and 40PTM; iRT and intensity) can place.
/// </summary>
/// <remarks>
/// mzLib declares these six models accept-all for UNIMOD ids: Koina's notes.yaml for the 2025 models says "All unimod
/// modifications are supported", and PTMs_gl runs the same preprocessing. The server still refuses some (residue, id)
/// pairs, and one refused peptide fails the whole request with HTTP 400. This class holds Koina's rule so mzLib can
/// refuse those peptides before sending.
///
/// Source: wilhelm-lab/koina at d17616997e45f6684a41a898c598e85d06130397. All six ensembles run the same two
/// preprocessing steps, Prosit_Preprocess_ac_gain and Prosit_Preprocess_ac_loss, so they share one rule:
///   https://github.com/wilhelm-lab/koina/blob/d17616997e45f6684a41a898c598e85d06130397/models/Prosit/Prosit_2025_irt_40PTM/config.pbtxt
///   (and Prosit_2025_irt_22PTM, Prosit_2025_intensity_22PTM, Prosit_2025_intensity_40PTM, Prosit_2024_irt_PTMs_gl and
///   Prosit_2024_intensity_PTMs_gl; at this sha the iRT 40PTM ensemble even calls the 22PTM core).
///
/// Koina encodes each modified position by the atoms it gains and loses. For each position it does, in order:
///   0. Looks the id up in its UNIMOD table, models/Deeplc/Deeplc_Preprocess_AC/1/unimod.obo (dated 2022-10-20,
///      newest id 2072). An id it lacks fails anywhere ("KeyError: 'UNIMOD:2073'").
///      Prosit_Preprocess_ac_gain/1/modifications.py, Unimod.lookup_sequence_m.
///   1. Uses the hand-written table for a trained pair: dict_ptm_atom_count_gain and dict_ptm_atom_count_loss in
///      models/Prosit/Prosit_Preprocess_ac_gain/1/model.py and models/Prosit/Prosit_Preprocess_ac_loss/1/model.py.
///      Both tables list the same 64 pairs. They are <see cref="TrainedTokens"/>.
///   2. Otherwise encodes the pair "zero-shot", from models/Prosit/Prosit_Preprocess_ac_gain/1/zero_shot.py
///      (identical in ac_loss): loss = the residue's attachment group, gain = loss + the UNIMOD delta composition.
///      It refuses the pair when the residue has no attachment group (<see cref="AttachmentGroups"/>), when the
///      composition holds an element other than H C N O P S or a symbol it cannot read, or when gain has a negative
///      atom count.
/// Koina's parser gives both termini the empty residue (modifications.py, ProformaParser.parse_sequence), so an
/// N-terminal and a C-terminal modification with the same id share one table row ("_1") and one attachment group.
/// That is why C-terminal tokens are accepted here: "PEPTIDEK-[UNIMOD:34]" and "-[UNIMOD:1]" return 200 live.
///
/// Tokens are written as KoinaModelBase writes them: "E[UNIMOD:954]" on a residue, "[UNIMOD:1]-" on the N-terminus
/// and "-[UNIMOD:2]" on the C-terminus.
///
/// The data below was generated from those files, not typed in, by a script that is not in this repository. Refresh it
/// whenever Koina changes Prosit_Preprocess_ac_gain/1/model.py, Prosit_Preprocess_ac_loss/1/model.py,
/// Prosit_Preprocess_ac_gain/1/zero_shot.py or models/Deeplc/Deeplc_Preprocess_AC/1/unimod.obo: diff those four files
/// against the sha above and regenerate. The script read the two tables, dict_attachment_loss and parse_composition's
/// constants with Python's ast.literal_eval (it never ran Koina's code), applied the same parsing to each
/// delta_composition in unimod.obo, and compared the result with mzLib's Omics/Resources/unimod.xml.
/// A check against every (position, id) pair for mzLib's 1,560 UNIMOD ids, 34,408 tokens, agreed with this class on all.
/// </remarks>
internal static class PrositAtomCountModificationRules
{
    /// <summary>
    /// The newest id in Koina's unimod.obo (2022-10-20). mzLib's unimod.xml has 47 newer ids (2073-2147) the server
    /// cannot look up. Every id up to 2072 in mzLib's catalog is also in Koina's file, so "known to Koina" is
    /// "in mzLib's catalog and at most 2072".
    /// </summary>
    internal const int NewestUnimodIdOnKoina = 2072;

    /// <summary>
    /// The trained pairs, from the keys of dict_ptm_atom_count_gain (the same keys as dict_ptm_atom_count_loss).
    /// A terminal key "_N" gives both "[UNIMOD:N]-" and "-[UNIMOD:N]". K_12118, K_19903 and K_129317 are left out,
    /// because those ids are not in unimod.obo and fail the lookup. K_5634 stays: lookup_sequence_m skips the lookup
    /// for UNIMOD:5634.
    /// </summary>
    internal static readonly IReadOnlySet<string> TrainedTokens = new HashSet<string>
    {
        "C[UNIMOD:4]", "C[UNIMOD:34]", "C[UNIMOD:35]", "C[UNIMOD:312]", "C[UNIMOD:1263]", "C[UNIMOD:2062]", "D[UNIMOD:34]", "E[UNIMOD:27]",
        "E[UNIMOD:34]", "H[UNIMOD:21]", "H[UNIMOD:34]", "H[UNIMOD:35]", "I[UNIMOD:34]", "K[UNIMOD:1]", "K[UNIMOD:4]", "K[UNIMOD:34]",
        "K[UNIMOD:35]", "K[UNIMOD:36]", "K[UNIMOD:37]", "K[UNIMOD:56]", "K[UNIMOD:58]", "K[UNIMOD:59]", "K[UNIMOD:121]", "K[UNIMOD:214]",
        "K[UNIMOD:535]", "K[UNIMOD:730]", "K[UNIMOD:737]", "K[UNIMOD:1263]", "K[UNIMOD:1289]", "K[UNIMOD:1293]", "K[UNIMOD:1848]", "K[UNIMOD:1990]",
        "K[UNIMOD:2016]", "K[UNIMOD:5634]", "L[UNIMOD:34]", "M[UNIMOD:35]", "N[UNIMOD:7]", "N[UNIMOD:34]", "P[UNIMOD:21]", "P[UNIMOD:35]",
        "Q[UNIMOD:7]", "Q[UNIMOD:28]", "Q[UNIMOD:34]", "R[UNIMOD:7]", "R[UNIMOD:34]", "R[UNIMOD:36]", "R[UNIMOD:267]", "S[UNIMOD:21]",
        "S[UNIMOD:43]", "T[UNIMOD:21]", "T[UNIMOD:43]", "W[UNIMOD:35]", "Y[UNIMOD:21]", "[UNIMOD:1]-", "[UNIMOD:58]-", "[UNIMOD:59]-",
        "[UNIMOD:214]-", "[UNIMOD:411]-", "[UNIMOD:730]-", "[UNIMOD:737]-", "[UNIMOD:2016]-", "-[UNIMOD:1]", "-[UNIMOD:58]", "-[UNIMOD:59]",
        "-[UNIMOD:214]", "-[UNIMOD:411]", "-[UNIMOD:730]", "-[UNIMOD:737]", "-[UNIMOD:2016]",
    };

    /// <summary>
    /// dict_attachment_loss in zero_shot.py: the atoms a residue loses when a modification attaches, keyed by residue,
    /// with "" for either terminus. A, E, F, G, N, P and V have none, so only their trained pairs are accepted.
    /// </summary>
    internal static readonly IReadOnlyDictionary<string, ChemicalFormula> AttachmentGroups = new Dictionary<string, string>
    {
        [""] = "H2N",
        ["C"] = "HS",
        ["D"] = "HO2",
        ["H"] = "HN",
        ["I"] = "H2N",
        ["K"] = "H2N",
        ["L"] = "H2N",
        ["M"] = "H3CS",
        ["Q"] = "H2N",
        ["R"] = "H2N",
        ["S"] = "HO",
        ["T"] = "HO",
        ["W"] = "H",
        ["Y"] = "HO",
    }.ToDictionary(group => group.Key, group => ChemicalFormula.ParseFormula(group.Value));

    /// <summary>
    /// Ids whose delta composition in unimod.obo zero_shot.parse_composition cannot read, although mzLib's formula
    /// for them holds only H C N O P S, so the element check below cannot see the problem. They are written with
    /// dHex (its token pattern needs an upper-case first letter), Me, Sulf, Ac, HexN or Hep. Most are fucosylated
    /// glycans. Koina refuses them on every residue it has no trained pair for.
    /// </summary>
    internal static readonly IReadOnlySet<int> UnreadableCompositionIds = new HashSet<int>
    {
        142, 145, 146, 147, 150, 152, 155, 156, 158, 295, 305, 307, 308, 454, 490, 910,
        1367, 1375, 1376, 1377, 1378, 1379, 1410, 1411, 1436, 1437, 1440, 1445, 1446, 1449, 1453, 1454,
        1459, 1462, 1463, 1467, 1471, 1472, 1474, 1475, 1476, 1477, 1478, 1481, 1484, 1485, 1491, 1492,
        1493, 1494, 1497, 1498, 1499, 1500, 1501, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1515,
        1516, 1518, 1519, 1520, 1522, 1523, 1525, 1526, 1527, 1528, 1529, 1533, 1534, 1535, 1536, 1537,
        1538, 1539, 1541, 1542, 1543, 1544, 1546, 1547, 1548, 1550, 1552, 1553, 1554, 1555, 1556, 1557,
        1562, 1564, 1571, 1579, 1580, 1581, 1584, 1586, 1587, 1588, 1593, 1594, 1596, 1597, 1606, 1607,
        1608, 1609, 1612, 1614, 1615, 1616, 1618, 1620, 1621, 1622, 1624, 1625, 1630, 1634, 1635, 1637,
        1639, 1640, 1641, 1643, 1644, 1645, 1648, 1650, 1651, 1652, 1653, 1654, 1656, 1657, 1658, 1660,
        1661, 1663, 1666, 1667, 1668, 1669, 1670, 1671, 1673, 1674, 1675, 1676, 1679, 1680, 1682, 1683,
        1684, 1686, 1688, 1689, 1691, 1693, 1695, 1697, 1698, 1699, 1701, 1702, 1703, 1706, 1707, 1708,
        1709, 1714, 1715, 1716, 1718, 1719, 1721, 1722, 1724, 1726, 1727, 1728, 1730, 1733, 1735, 1736,
        1737, 1739, 1740, 1743, 1746, 1747, 1748, 1750, 1751, 1754, 1758, 1760, 1761, 1762, 1764, 1765,
        1766, 1767, 1768, 1770, 1771, 1774, 1775, 1778, 1781, 1782, 1783, 1784, 1785, 1786, 1840, 1933,
        1935, 1937, 1941, 1943, 1946, 1948, 1949, 1950, 1951, 1952, 1954, 1955, 1956, 1958, 1960, 1963,
        1965, 1966, 1967, 1969,
    };

    private static readonly string[] EncodedElements = { "H", "C", "N", "O", "P", "S" };
    private static readonly Regex ElementSymbol = new(@"[A-Z][a-z]?", RegexOptions.Compiled);
    private static readonly Regex TokenPattern = new(
        @"^(?:(?<residue>[A-Z])\[UNIMOD:(?<id>\d+)\]|\[UNIMOD:(?<id>\d+)\]-|-\[UNIMOD:(?<id>\d+)\])$", RegexOptions.Compiled);

    /// <summary>
    /// Each UNIMOD id's delta composition, from mzLib's own catalog. For every id Koina can read, these formulas
    /// (isotopes counted as their element) give the same H C N O P S counts as unimod.obo, so the checks below see
    /// what the server sees.
    /// </summary>
    private static readonly Lazy<IReadOnlyDictionary<int, ChemicalFormula>> UnimodFormulas = new(() =>
        Mods.UnimodModifications
            .Where(mod => mod.ChemicalFormula != null && mod.DatabaseReference != null
                && mod.DatabaseReference.TryGetValue("Unimod", out var ids) && ids.Count > 0 && int.TryParse(ids[0], out _))
            .GroupBy(mod => int.Parse(mod.DatabaseReference["Unimod"][0]))
            .ToDictionary(group => group.Key, group => group.First().ChemicalFormula));

    /// <summary>
    /// Why Koina refuses this modification token, or null when it accepts it.
    /// </summary>
    /// <param name="token">"E[UNIMOD:954]" on a residue, "[UNIMOD:1]-" on the N-terminus, "-[UNIMOD:2]" on the C-terminus.</param>
    internal static string? RefusalReason(string token)
    {
        var match = TokenPattern.Match(token);
        if (!match.Success)
            throw new ArgumentException($"'{token}' is not a modification token like \"E[UNIMOD:954]\", \"[UNIMOD:1]-\" or \"-[UNIMOD:2]\".", nameof(token));

        var position = match.Groups["residue"].Success ? ModificationPositionType.Residue
            : token[0] == '[' ? ModificationPositionType.NTerminus : ModificationPositionType.CTerminus;
        return RefusalReason(position, match.Groups["residue"].Value, int.Parse(match.Groups["id"].Value));
    }

    /// <summary>
    /// Why Koina refuses this modification, where it will be sent, or null when it accepts it.
    /// </summary>
    internal static string? Decide(CanonicalSequence sequence, CanonicalModification modification)
    {
        string residue = modification.PositionType == ModificationPositionType.Residue
            ? sequence.BaseSequence[modification.ResidueIndex!.Value].ToString()
            : "";
        return RefusalReason(modification.PositionType, residue, modification.UnimodId!.Value);
    }

    /// <summary>
    /// The rule for one position; <paramref name="residue"/> is "" on either terminus, as Koina keys it.
    /// </summary>
    private static string? RefusalReason(ModificationPositionType position, string residue, int id)
    {
        string token = position switch
        {
            ModificationPositionType.NTerminus => $"[UNIMOD:{id}]-",
            ModificationPositionType.CTerminus => $"-[UNIMOD:{id}]",
            _ => $"{residue}[UNIMOD:{id}]"
        };
        if (TrainedTokens.Contains(token))
            return null;

        if (id > NewestUnimodIdOnKoina || !UnimodFormulas.Value.TryGetValue(id, out var delta))
            return $"UNIMOD:{id} is not in Koina's UNIMOD table";

        if (!AttachmentGroups.TryGetValue(residue, out var group))
            return $"residue {residue} has no attachment group, and {token} is not a trained modification";

        var otherElements = ElementSymbol.Matches(delta.Formula).Select(symbol => symbol.Value)
            .Where(symbol => !EncodedElements.Contains(symbol)).Distinct().ToList();
        if (otherElements.Count > 0)
            return $"UNIMOD:{id} contains {string.Join(", ", otherElements)}, and Koina encodes only H, C, N, O, P and S";

        if (UnreadableCompositionIds.Contains(id))
            return $"Koina cannot read the composition of UNIMOD:{id} (written with dHex, Me, Sulf, Ac, HexN or Hep)";

        var negative = EncodedElements.Where(element => group.CountWithIsotopes(element) + delta.CountWithIsotopes(element) < 0).ToList();
        if (negative.Count > 0)
        {
            string site = position switch
            {
                ModificationPositionType.NTerminus => "the N-terminus",
                ModificationPositionType.CTerminus => "the C-terminus",
                _ => $"residue {residue}"
            };
            return $"the attachment group of {site} ({group.Formula}) plus UNIMOD:{id} ({delta.Formula}) has a negative {string.Join(", ", negative)} count";
        }

        return null;
    }
}
