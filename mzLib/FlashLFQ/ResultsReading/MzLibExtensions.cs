using Readers.ExternalResults.BaseClasses;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public static class MzLibExtensions
    {
        /// <summary>
        /// Makes a list of identification objects usable by FlashLFQ from an IQuantifiableResultFile
        /// </summary>
        public static List<Identification> MakeIdentifications(this IQuantifiableResultFile quantifiable)
        {
            IEnumerable<IQuantifiableRecord> quantifiableRecords = quantifiable.GetQuantifiableResults();
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
                if (allFiles.TryGetValue(record.FileName, out var fileInfo))
                {
                    // placeholder values for SpectraFileInfo that will be edited later
                    file = new SpectraFileInfo(record.FileName, "", 1, 1, 1);
                }
                else
                {
                    file = new SpectraFileInfo(record.FileName, "", 1, 1, 1);
                    allFiles.Add(record.FileName, fileInfo);
                }

                List<ProteinGroup> proteinGroups = new();
                foreach (var info in record.ProteinGroupInfos)
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