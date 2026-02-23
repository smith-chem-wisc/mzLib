using System;
using System.Collections.Generic;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;

namespace Development.QuantificationDevelopment.TestHelpers;

/// <summary>
/// Hardcoded experimental design for the SpikeIn-5mix-MS3 TMT dataset.
/// 14 files × 10 channels (TMT10-plex) = 140 IsobaricQuantSampleInfo objects.
/// UPS1 spike-in concentrations (0.125, 0.5, 0.667, 1.0) vary by mixture and channel.
/// Channels 126 and 131 are pooled reference channels in every file.
/// </summary>
public class SpikeInExperimentalDesign : IExperimentalDesign
{
    public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

    public SpikeInExperimentalDesign()
    {
        FileNameSampleInfoDictionary = BuildDesign();
    }

    private static Dictionary<string, ISampleInfo[]> BuildDesign()
    {
        var design = new Dictionary<string, ISampleInfo[]>();

        // TMT10-plex channel labels and reporter ion m/z values (in order)
        string[] channelLabels = { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N" };
        double[] channelMzs   = { 126.12776, 127.12476, 127.13108, 128.12811, 128.13443, 129.13147, 129.13779, 130.13482, 130.14114, 131.13818 };

        // Channel-to-condition mappings per mixture type.
        (string condition, int bioRep)[] mixture15Conditions =
        {
            ("Reference", 1), // 126
            ("0.667",     1), // 127N
            ("0.125",     1), // 127C
            ("0.5",       1), // 128N
            ("1",         1), // 128C
            ("0.125",     2), // 129N
            ("0.5",       2), // 129C
            ("1",         2), // 130N
            ("0.667",     2), // 130C
            ("Reference", 1), // 131
        };

        (string condition, int bioRep)[] mixture2Conditions =
        {
            ("Reference", 1), // 126
            ("0.5",       1), // 127N
            ("1",         1), // 127C
            ("0.667",     1), // 128N
            ("0.125",     1), // 128C
            ("1",         2), // 129N
            ("0.667",     2), // 129C
            ("0.125",     2), // 130N
            ("0.5",       2), // 130C
            ("Reference", 1), // 131
        };

        (string condition, int bioRep)[] mixture3Conditions =
        {
            ("Reference", 1), // 126
            ("0.125",     1), // 127N
            ("0.667",     1), // 127C
            ("1",         1), // 128N
            ("0.5",       1), // 128C
            ("0.5",       2), // 129N
            ("0.125",     2), // 129C
            ("0.667",     2), // 130N
            ("1",         2), // 130C
            ("Reference", 1), // 131
        };

        (string condition, int bioRep)[] mixture4Conditions =
        {
            ("Reference", 1), // 126
            ("1",         1), // 127N
            ("0.5",       1), // 127C
            ("0.125",     1), // 128N
            ("0.667",     1), // 128C
            ("0.667",     2), // 129N
            ("1",         2), // 129C
            ("0.5",       2), // 130N
            ("0.125",     2), // 130C
            ("Reference", 1), // 131
        };

        var files = new[]
        {
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01", 1, 1),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02", 1, 2),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03", 1, 3),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_01", 2, 1),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_02", 2, 2),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_03", 2, 3),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_01", 3, 1),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_02", 3, 2),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_03", 3, 3),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_02", 4, 2),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_03", 4, 3),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_01", 5, 1),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_02", 5, 2),
            ("161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_03", 5, 3),
        };

        foreach (var (fileName, mixtureNum, techRep) in files)
        {
            var conditions = mixtureNum switch
            {
                1 or 5 => mixture15Conditions,
                2      => mixture2Conditions,
                3      => mixture3Conditions,
                4      => mixture4Conditions,
                _      => throw new InvalidOperationException($"Unknown mixture number: {mixtureNum}")
            };

            design[fileName] = CreateChannelInfos(fileName, mixtureNum, techRep, channelLabels, channelMzs, conditions);
        }

        return design;
    }

    private static ISampleInfo[] CreateChannelInfos(
        string fileName,
        int mixtureNumber,
        int techRep,
        string[] channelLabels,
        double[] channelMzs,
        (string condition, int bioRep)[] channelConditions)
    {
        var infos = new ISampleInfo[channelLabels.Length];
        for (int i = 0; i < channelLabels.Length; i++)
        {
            bool isRef = channelConditions[i].condition == "Reference";
            infos[i] = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: fileName,
                condition: channelConditions[i].condition,
                biologicalReplicate: channelConditions[i].bioRep,
                technicalReplicate: techRep,
                fraction: 0,
                plexId: mixtureNumber,
                channelLabel: channelLabels[i],
                reporterIonMz: channelMzs[i],
                isReferenceChannel: isRef);
        }
        return infos;
    }
}
