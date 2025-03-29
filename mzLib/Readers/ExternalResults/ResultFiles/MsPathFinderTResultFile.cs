using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;

namespace Readers
{
    public class MsPathFinderTResultFile : ResultFile<MsPathFinderTResult>, IResultFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public MsPathFinderTResultFile(string filePath) : base(filePath, Software.MsPathFinderT) //instantiates w/ SupportedFileType ending
        {
            FileType = filePath.ParseFileType();
        }

        public MsPathFinderTResultFile() : base()
        {
            FileType = FilePath.IsNullOrEmpty() ? SupportedFileType.MsPathFinderTAllResults : FilePath.ParseFileType();
        }

        public override string ToString() //TODO consider removing
        {
            StringBuilder x = new StringBuilder();
            
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();
            foreach (MsPathFinderTResult res in Results)
            {
                x.Append(res.BaseSequence);
            }

            
            return x.ToString();
        }

        public Dictionary<string, string> ToDict() //key,value -> base,record pair TODO consider removing.
            // Here's what a single value entry looks like: {[SLEVFEKLEAKVQQAIDTITLLQMEIEELKEKNNSLSQEVQNAQHQREELERENNHLKEQQNGWQERLQALLGRMEEV, 1955	M	SLEVFEKLEAKVQQAIDTITLLQMEIEELKEKNNSLSQEVQNAQHQREELERENNHLKEQQNGWQERLQALLGRMEEV	-		C396H647N117O135S2 : 9265.680070361872	sp|P0AF36|ZAPB_ECOLI	Cell division protein ZapB OS=Escherichia coli (strain K12) OX=83333 GN=zapB PE=1 SV=1	82	4	81	10	928.0769611	9265.680072	0	65	1	9.99E-308	9.99E-308	0	0	]}
        {
            var dict = new Dictionary<string, string>();
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();

            foreach (MsPathFinderTResult res in Results)
            {
                if (!string.IsNullOrEmpty(res.BaseSequence) && !dict.ContainsKey(res.BaseSequence))
                {
                    StringBuilder SB = new StringBuilder(); 
                    dict[res.BaseSequence] = 
                        SB.Append(res.OneBasedScanNumber)
                          .Append("\t") 
                           .Append(res.PreviousResidue)
                           .Append("\t")
                           .Append(res.BaseSequence)
                           .Append("\t")
                           .Append(res.NextResidue)
                           .Append("\t")
                           .Append(res.Modifications)
                           .Append("\t")
                           .Append(res.ChemicalFormula)
                           .Append("\t")
                           .Append(res.ProteinName)
                           .Append("\t")
                           .Append(res.ProteinDescription)
                           .Append("\t")
                           .Append(res.Length)
                           .Append("\t")
                           .Append(res.OneBasedStartResidue)
                           .Append("\t")
                           .Append(res.OneBasedEndResidue)
                           .Append("\t")
                           .Append(res.Charge)
                           .Append("\t")
                           .Append(res.MostAbundantIsotopeMz)
                           .Append("\t")
                           .Append(res.MonoisotopicMass)
                           .Append("\t")
                           .Append(res.Ms1Features)
                           .Append("\t")
                           .Append(res.NumberOfMatchedFragments)
                           .Append("\t")
                           .Append(res.Probability)
                           .Append("\t")
                           .Append(res.SpecEValue)
                           .Append("\t")
                           .Append(res.EValue)
                           .Append("\t")
                           .Append(res.QValue)
                           .Append("\t")
                           .Append(res.PepQValue)
                           .Append("\t")
                           .Append(res.FileNameWithoutExtension)
                           .ToString();
                }
            }

            return dict;
        }
        public Dictionary<string, List<string>> ToDictList() // TODO consider removing.
        {
            var dictList = new Dictionary<string, List<string>>();
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();

            foreach (MsPathFinderTResult res in Results)
            {
                string uniqueKey = $"{res.BaseSequence}_{res.Modifications}";
                if (!string.IsNullOrEmpty(uniqueKey) && !dictList.ContainsKey(uniqueKey))
                {
                    dictList[uniqueKey] = new List<string>
                    {
                        res.OneBasedScanNumber.ToString(),
                        res.PreviousResidue.ToString(),
                        res.BaseSequence,
                        res.NextResidue.ToString(),
                        res.Modifications,
                        res.ChemicalFormula.ToString(),
                        res.ProteinName,
                        res.ProteinDescription,
                        res.Length.ToString(),
                        res.OneBasedStartResidue.ToString(),
                        res.OneBasedEndResidue.ToString(),
                        res.Charge.ToString(),
                        res.MostAbundantIsotopeMz.ToString(),
                        res.MonoisotopicMass.ToString(),
                        res.Ms1Features.ToString(),
                        res.NumberOfMatchedFragments.ToString(),
                        res.Probability.ToString(),
                        res.SpecEValue.ToString(),
                        res.EValue.ToString(),
                        res.QValue.ToString(),
                        res.PepQValue.ToString(),
                        res.FileNameWithoutExtension
                    };
                }
            }

            return dictList;
        }

        public Dictionary<string, List<List<string>>> ToDictListList()
        {
            var dictList = new Dictionary<string, List<List<string>>>();
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();

            foreach (MsPathFinderTResult res in Results)
            {
                string baseKey = res.BaseSequence;
                var valuesList = new List<string>
                {
                    res.OneBasedScanNumber.ToString(),
                    res.PreviousResidue.ToString(),
                    res.BaseSequence,
                    res.NextResidue.ToString(),
                    res.Modifications,
                    res.ChemicalFormula.ToString(),
                    res.ProteinName,
                    res.ProteinDescription,
                    res.Length.ToString(),
                    res.OneBasedStartResidue.ToString(),
                    res.OneBasedEndResidue.ToString(),
                    res.Charge.ToString(),
                    res.MostAbundantIsotopeMz.ToString(),
                    res.MonoisotopicMass.ToString(),
                    res.Ms1Features.ToString(),
                    res.NumberOfMatchedFragments.ToString(),
                    res.Probability.ToString(),
                    res.SpecEValue.ToString(),
                    res.EValue.ToString(),
                    res.QValue.ToString(),
                    res.PepQValue.ToString(),
                    res.FileNameWithoutExtension
                };

                if (!dictList.ContainsKey(baseKey))
                {
                    dictList[baseKey] = new List<List<string>>();
                }

                dictList[baseKey].Add(valuesList);
            }

            return dictList;
        }


        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsPathFinderTResult.CsvConfiguration);
            Results = csv.GetRecords<MsPathFinderTResult>().ToList();
            if (Results.Any() && Results.First().FileNameWithoutExtension.IsNullOrEmpty())
                Results.ForEach(p => p.FileNameWithoutExtension = string.Join("_", Path.GetFileNameWithoutExtension(FilePath).Split('_')[..^1]));
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsPathFinderTResult.CsvConfiguration))
            {
                csv.WriteHeader<MsPathFinderTResult>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }
}
