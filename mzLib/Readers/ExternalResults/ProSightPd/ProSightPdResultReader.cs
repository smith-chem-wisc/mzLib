using System.Data.SQLite;

namespace Readers;

public static class ProSightPdResultReader
{
    public static List<ProSightPdPsmRecord> LoadPsms(string sqlitePath)
    {
        using var connection = new SQLiteConnection($"Data Source={sqlitePath};Version=3;");
        connection.Open();

        using var command = new SQLiteCommand("SELECT * FROM TargetPsms;", connection);
        using var reader = command.ExecuteReader();

        var results = new List<ProSightPdPsmRecord>();
        while (reader.Read())
        {
            results.Add(new ProSightPdPsmRecord
            {
                WorkflowID = reader.GetInt32(0),
                ID = reader.GetInt32(1),
                ActivationTypes = reader.GetInt32(2),
                IonInjectTime = reader.IsDBNull(3) ? null : reader.GetDouble(3),
                Intensity = reader.GetDouble(4),
                RetentionTime = reader.GetDouble(5),
                FragmentationScans = reader.IsDBNull(6) ? null : reader.GetString(6),
                NumberFragmentationScans = reader.GetInt32(7),
                PrecursorScans = reader.IsDBNull(8) ? null : reader.GetString(8),
                NumberPrecursorScans = reader.GetInt32(9),
                MSOrderTypes = reader.GetInt32(10),
                IdentifyingNodeTypeName = reader.IsDBNull(11) ? null : reader.GetString(11),
                IdentifyingNodeName = reader.IsDBNull(12) ? null : reader.GetString(12),
                IdentifyingNodeNumber = reader.GetInt32(13),
                IdentifyingNodeSearchID = reader.IsDBNull(14) ? null : reader.GetString(14),
                UniqueSequenceID = reader.GetInt32(15),
                Sequence = reader.IsDBNull(16) ? null : reader.GetString(16),
                ModifiedSequence = reader.IsDBNull(17) ? null : reader.GetString(17),
                UniqueSequence = reader.IsDBNull(18) ? null : reader.GetString(18),
                MatchConfidence = reader.GetInt32(19),
                Modifications = reader.IsDBNull(20) ? null : reader.GetString(20),
                ParentProteinCount = reader.GetInt32(21),
                ParentProteinAccessions = reader.IsDBNull(22) ? null : reader.GetString(22),
                ParentProteinDescriptions = reader.IsDBNull(23) ? null : reader.GetString(23),
                Rank = reader.GetInt32(24),
                SearchEngineRank = reader.GetInt32(25),
                Charge = reader.GetInt32(26),
                OriginalPrecursorCharge = reader.GetInt32(27),
                MassOverCharge = reader.GetDouble(28),
                DetectedNeutralMass = reader.GetDouble(29),
                TheoreticalNeutralMass = reader.GetDouble(30),
                DeltaMassInDa = reader.GetDouble(31),
                DeltaMassInPpm = reader.GetDouble(32),
                IonsMatched = reader.IsDBNull(33) ? null : reader.GetValue(33).ToString(),
                MatchedIonsCount = reader.GetInt32(34),
                TotalIonsCount = reader.GetInt32(35),
                SpectrumFileId = reader.GetInt32(36),
                StudyFileId = reader.IsDBNull(37) ? null : reader.GetString(37),
                ExcludedBy = reader.GetInt32(38),
                CScore = reader.IsDBNull(39) ? null : reader.GetDouble(39),
                ResidueCleavages = reader.IsDBNull(40) ? null : reader.GetDouble(40),
                CorrectedDeltaMassDa = reader.IsDBNull(41) ? null : reader.GetDouble(41),
                CorrectedDeltaMassppm = reader.IsDBNull(42) ? null : reader.GetDouble(42),
                DetectedIonCount = reader.GetInt32(43),
                ProteoformLevel = reader.GetInt32(44),
                ExternalTopDownDisplays = reader.IsDBNull(45) ? null : reader.GetString(45),
                CompensationVoltage = reader.IsDBNull(46) ? null : reader.GetDouble(46),
                PTMsLocalized = reader.GetInt32(47),
                PTMsIdentified = reader.GetInt32(48),
                SequenceDefined = reader.GetInt32(49),
                GeneIdentified = reader.GetInt32(50),
                ProformaHash = reader.IsDBNull(51) ? null : reader.GetString(51),
                LocalizedModifications = reader.IsDBNull(52) ? null : reader.GetString(52),
                NonlocalizedModifications = reader.IsDBNull(53) ? null : reader.GetString(53),
                NCE = reader.IsDBNull(54) ? null : reader.GetDouble(54),
                LogPScore = reader.IsDBNull(55) ? null : reader.GetDouble(55),
                LogEValue = reader.IsDBNull(56) ? null : reader.GetDouble(56),
                SubSeqAARange = reader.IsDBNull(57) ? null : reader.GetString(57),
                Qvalue = reader.GetDouble(58),
            });
        }

        return results;
    }
}
