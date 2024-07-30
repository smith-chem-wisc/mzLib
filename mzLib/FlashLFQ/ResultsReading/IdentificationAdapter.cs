using Readers.ExternalResults.BaseClasses;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.ResultsReading
{
    public static class IdentificationAdapter
    {
        public static List<Identification> MakeIdentifications(IEnumerable<IQuantifiableRecord> quantifiableRecords)
        {
            List<Identification> identifications = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            Dictionary<string, SpectraFileInfo> allFiles = new Dictionary<string, SpectraFileInfo>();

            foreach (var record in quantifiableRecords)
            {
                string baseSequence = record.BaseSequence;
                string modifiedSequence = record.ModifiedSequence;
                double ms2RetentionTimeInMinutes = record.RetentionTime;
                double monoisotopicMass = record.MonoisotopicMass;
                int precursurChargeState = record.ChargeState;

                SpectraFileInfo file = null;
                if (allFiles.TryGetValue(record.FullFilePath, out var fileInfo))
                {
                    file = new SpectraFileInfo(record.FullFilePath, "", 1, 1, 1);
                }
                else
                {
                    file = new SpectraFileInfo(record.FullFilePath, "", 1, 1, 1);
                    allFiles.Add(record.FullFilePath, fileInfo);
                }

                List<ProteinGroup> proteinGroups = new();
                foreach (var info in record.proteinGroupInfos)
                {
                    if (allProteinGroups.TryGetValue(info.proteinAccessions, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(info.proteinAccessions, new ProteinGroup(info.proteinAccessions, info.geneName, info.organism));
                        proteinGroups.Add(allProteinGroups[info.proteinAccessions]);
                    }
                }
                Identification id = new Identification(file, baseSequence, modifiedSequence, monoisotopicMass, ms2RetentionTimeInMinutes, precursurChargeState, proteinGroups);
                identifications.Add(id);

            }

            return identifications;
        }
    }
}
