using MzLibUtil;
using PredictionClients.Koina.AbstractClasses;
using System.ComponentModel;
using System.Text;

namespace PredictionClients.Koina.SupportedModels.RetentionTimeModels
{
    public class Prosit2020iRTTMT : RetentionTimeModel
    {
        public override string ModelName => "Prosit_2020_irt_TMT";
        public override int MaxBatchSize => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 1;
        public override bool IsIndexedRetentionTimeModel => true;
        public override Dictionary<string, string> ValidModificationUnimodMapping => new()
        {
            {"[Common Variable:Oxidation on M]", "[UNIMOD:35]"},
            {"[Common Fixed:Carbamidomethyl on C]", "[UNIMOD:4]"},
            // SILAC modifications
            {"[Common Variable:Label:13C(6)15N(2) on K]", "[UNIMOD:259]"},
            {"[Common Variable:Label:13C(6)15N(4) on R]", "[UNIMOD:267]"},
            // TMT modifications
            {"[Common Fixed:TMT6plex on K]", "[UNIMOD:737]"},
            {"[Common Fixed:TMT6plex on N-terminus]", "[UNIMOD:737]-"},
            // TMTpro modifications
            {"[Common Fixed:TMTpro on K]", "[UNIMOD:2016]"},
            {"[Common Fixed:TMTpro on N-terminus]", "[UNIMOD:2016]-"},
            // iTRAQ 4-plex modifications
            {"[Common Fixed:iTRAQ4plex on K]", "[UNIMOD:214]"},
            {"[Common Fixed:iTRAQ4plex on N-terminus]", "[UNIMOD:214]-"},
            // iTRAQ 8-plex modifications
            {"[Common Fixed:iTRAQ8plex on K]", "[UNIMOD:730]"},
            {"[Common Fixed:iTRAQ8plex on N-terminus]", "[UNIMOD:730]-"}
        };
        public override List<string> PeptideSequences { get; } = new();
        public override List<PeptideRTPrediction> Predictions { get; protected set; } = new();


        public Prosit2020iRTTMT(List<string> peptideSequences, out WarningException? warnings)
        {
            if (peptideSequences.IsNullOrEmpty())
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            //Ensure all peptides are valid
            var invalidSequences = new List<string>();
            foreach (var seq in peptideSequences)
            {
                if (!IsValidSequence(seq) || !HasValidModifications(seq))
                {
                    invalidSequences.Add(seq);
                }
                else
                {
                    PeptideSequences.Add(ConvertToPrositModificationFormat(seq));
                }
            }

            warnings = null;
            if (invalidSequences.Count > 0)
            {
                var sb = new StringBuilder();
                sb.AppendLine("The following peptide sequences were invalid and will be skipped:");
                foreach (var invalid in invalidSequences)
                {
                    sb.AppendLine(invalid);
                }
                warnings = new WarningException(sb.ToString());
            }
        }
        public override List<Dictionary<string, object>> ToBatchedRequests()
        {
            var batchedPeptides = PeptideSequences.Chunk(MaxBatchSize).ToList();
            var batchedRequests = new List<Dictionary<string, object>>();
            for (int i = 0; i < batchedPeptides.Count; i++)
            {
                var request = new Dictionary<string, object>
                {
                    { "id", $"Batch{i}_" + Guid.NewGuid()},
                    { "inputs", new List<object>
                        {
                            new {
                                name = "peptide_sequences",
                                shape = new[]{ batchedPeptides[i].Length, 1 },
                                datatype = "BYTES",
                                data = batchedPeptides[i]
                            },
                        }
                    }
                };
                batchedRequests.Add(request);
            }
            return batchedRequests;
        }
    }
}
