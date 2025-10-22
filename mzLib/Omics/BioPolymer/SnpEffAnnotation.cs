using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Specifications are described here: http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
    /// </summary>
    public class SnpEffAnnotation
    {
        private static readonly Regex HGVSProteinRegex = new Regex(@"(p\.)([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])");

        // All public getters: ensure they are always initialized (never left unassigned).
        /// <summary>
        /// Original SnpEff annotation string.
        /// </summary>
        public string Annotation { get; }
        public string Allele { get; } = string.Empty;
        public string[] Effects { get; } = Array.Empty<string>();
        public string PutativeImpact { get; } = string.Empty;
        public string GeneName { get; } = string.Empty;
        public string GeneID { get; } = string.Empty;

        /// <summary>
        /// It looks like these are sometimes domains, like the ones annotated in UniProt,
        /// Otherwise, this tends to just be "transcript"
        ///
        /// Some examples:
        /// sequence_feature: can be initiator-methionine:Removed ... maybe not too helpful for proteomics, since this is assumed
        /// sequence_feature: helix:combinatorial_evidence_used_in_manual_assertion
        /// sequence_feature: nucleotide-phosphate-binding-region:ATP
        /// sequence_feature: domain:EGF-like_2
        /// sequence_feature: transmembrane-region:Transmembrane_region
        /// sequence_feature: topological-domain:Extracellular
        /// sequence_feature: modified-residue:phosphoserine
        /// </summary>
        public string FeatureType { get; } = string.Empty;
        /// <summary>
        /// Always seems to be the transcriptID
        /// </summary>
        public string FeatureID { get; } = string.Empty;
        public string TranscriptBiotype { get; } = string.Empty;
        public int ExonIntronRank { get; }
        public int ExonIntronTotal { get; }
        public string HGVSNotationDnaLevel { get; } = string.Empty;// kind of bad for ins and del because they notation aligns to most 3' coordinate, rather than leftmost
        public string HGVSNotationProteinLevel { get; } = string.Empty;
        public int OneBasedTranscriptCDNAPosition { get; }
        public int TranscriptCDNALength { get; }
        public int OneBasedCodingDomainSequencePosition { get; }
        public int CodingDomainSequenceLengthIncludingStopCodon { get; }
        public int OneBasedProteinPosition { get; }
        public int ProteinLength { get; }
        /// <summary>
        /// up/downstream: distance to first / last codon
        /// intergenic: distance to closest gene
        /// exonic: distance to closest intron boundary (+ is upstream, - is downstream)
        /// intronic: distance to closest exon boundary (+ is upstream, - is downstream)
        /// motif: distance to first base in MOTIF
        /// miRNA: distance to first base in miRNA
        /// splice_site: distance to exon-intron boundary
        /// splice_region: distance to exon-intron boundary
        /// chip seq peak: distance to summit or peak center
        /// histone mark/state: distance to summit or peak center
        /// </summary>
        public int DistanceToFeature { get; }
        public string[] Warnings { get; } = Array.Empty<string>();

        public int AminoAcidLocation { get; }
        public char ReferenceAminoAcid { get; }
        public char AlternateAminoAcid { get; }
        public bool Missense { get; }
        public bool Synonymous { get; }
        public bool FrameshiftVariant { get; }
        public bool BadTranscript { get; }

        public SnpEffAnnotation(string annotation)
        {
            bool isSnpEffAnnotation = annotation.StartsWith("ANN=") || annotation.StartsWith("EFF=");
            Annotation = isSnpEffAnnotation ? annotation.Substring(4) : annotation;

            // If not a recognized snpEff style annotation, leave defaults (all properties already initialized)
            if (!isSnpEffAnnotation)
            {
                return;
            }

            // Split safely. Minimal examples (e.g. ANN=X|Y) produce few tokens.
            string[] a = Annotation.Split('|');

            string Get(int idx) => idx >= 0 && idx < a.Length ? a[idx] : string.Empty;

            Allele = Get(0);
            var effectsField = Get(1);
            Effects = string.IsNullOrEmpty(effectsField)
                ? Array.Empty<string>()
                : effectsField.Split('&', StringSplitOptions.RemoveEmptyEntries);

            PutativeImpact = Get(2);
            GeneName = Get(3);
            GeneID = Get(4);
            FeatureType = Get(5);
            FeatureID = Get(6);
            TranscriptBiotype = Get(7);

            // Exon/Intron rank/total: field 8 (e.g. "3/12")
            var exonIntron = Get(8);
            if (!string.IsNullOrEmpty(exonIntron))
            {
                var parts = exonIntron.Split('/');
                if (parts.Length > 0 && int.TryParse(parts[0], out int x)) ExonIntronRank = x;
                if (parts.Length > 1 && int.TryParse(parts[1], out int y)) ExonIntronTotal = y;
            }

            HGVSNotationDnaLevel = Get(9);
            HGVSNotationProteinLevel = Get(10);

            void ParseSlashField(string value, ref int first, ref int second)
            {
                if (string.IsNullOrEmpty(value)) return;
                var parts = value.Split('/');
                if (parts.Length > 0 && int.TryParse(parts[0], out int x)) first = x;
                if (parts.Length > 1 && int.TryParse(parts[1], out int y)) second = y;
            }

            {
                int pos = OneBasedTranscriptCDNAPosition;
                int len = TranscriptCDNALength;
                ParseSlashField(Get(11), ref pos, ref len);
                OneBasedTranscriptCDNAPosition = pos;
                TranscriptCDNALength = len;
            }
            {
                int pos = OneBasedCodingDomainSequencePosition;
                int len = CodingDomainSequenceLengthIncludingStopCodon;
                ParseSlashField(Get(12), ref pos, ref len);
                OneBasedCodingDomainSequencePosition = pos;
                CodingDomainSequenceLengthIncludingStopCodon = len;
            }
            {
                int pos = OneBasedProteinPosition;
                int len = ProteinLength;
                ParseSlashField(Get(13), ref pos, ref len);
                OneBasedProteinPosition = pos;
                ProteinLength = len;
            }

            if (int.TryParse(Get(14), out int dist))
            {
                DistanceToFeature = dist;
            }

            var warningsField = Get(15);
            Warnings = string.IsNullOrEmpty(warningsField)
                ? Array.Empty<string>()
                : warningsField.Split('&', StringSplitOptions.RemoveEmptyEntries);

            // Derive flags based on Effects (safe even if empty)
            Missense = Effects.Any(eff => eff == "missense_variant");
            FrameshiftVariant = Effects.Contains("frameshift_variant");

            Synonymous = Effects.Length == 0
                ? false // With no effect terms, treat as non-synonymous=false, synonymous=false (neutral/unknown)
                : !Effects.Any(eff => NonSynonymousVariations.Contains(eff));

            BadTranscript = Warnings.Any(w => BadTranscriptWarnings.Contains(w));

            // Additional amino acid / HGVS-level fields (if needed in future) can be derived here.
            // For now, keep defaults (0 / '\0').
        }

        //NOTE: The following arrays are retained for reference, but not currently used.

        //private string[] HighPutativeImpactEffects = new string[]
        //{
        //    "chromosome_number_variation",
        //    "exon_loss_variant",
        //    "frameshift_variant",
        //    "rare_amino_acid_variant",
        //    "splice_acceptor_variant", // often with intron_variant, sometimes with splice_donor_variant
        //    "splice_donor_variant", // often with intron_variant, sometimes with splice_acceptor_variant
        //    "start_lost",
        //    "stop_gained",
        //    "stop_lost",
        //    "transcript_ablation",
        //};

        //private string[] ModeratePutativeImpactEffects = new string[]
        //{
        //    "3_prime_UTR_truncation", "exon_loss", // appear together
        //    "5_prime_UTR_truncation", "exon_loss_variant", // appear together
        //    "coding_sequence_variant", // not seen much? Probably because missense is used more often.
        //    "conservative_inframe_insertion",
        //    "conservative_inframe_deletion",
        //    "disruptive_inframe_deletion",
        //    "disruptive_inframe_insertion",
        //    "inframe_deletion",// not common, in favor of more specific terms above
        //    "inframe_insertion",// not common, in favor of more specific terms above
        //    "missense_variant",
        //    "regulatory_region_ablation", // not common?
        //    "splice_region_variant", // often combined with intron_variant and non_coding_transcript_exon_variant
        //    "TFBS_ablation", // not common?
        //};

        private string[] NonSynonymousVariations = new string[]
        {
            "exon_loss_variant",
            "frameshift_variant",
            "rare_amino_acid_variant",
            "start_lost",
            "stop_gained",
            "stop_lost",
            "conservative_inframe_insertion",
            "conservative_inframe_deletion",
            "disruptive_inframe_deletion",
            "disruptive_inframe_insertion",
            "inframe_deletion", // not common, in favor of more specific terms above
            "inframe_insertion", // not common, in favor of more specific terms above
            "missense_variant",
        };

        //NOTE: The following arrays are retained for reference, but not currently used.

        //private string[] LowPutativeImpactEffects = new string[]
        //{
        //    "5_prime_UTR_premature_start_codon_gain_variant",
        //    "initiator_codon_variant",
        //    "splice_region_variant",
        //    "start_retained", // not used in human, with only one canonical start codon
        //    "stop_retained_variant", // fairly common
        //    "synonymous_variant",
        //    "sequence_feature"
        //};

        //private string[] ModifierEffects = new string[]
        //{
        //    "3_prime_UTR_variant",
        //    "5_prime_UTR_variant",
        //    "coding_sequence_variant",
        //    "conserved_intergenic_variant",
        //    "conserved_intron_variant",
        //    "downstream_gene_variant",
        //    "exon_variant",
        //    "feature_elongation",
        //    "feature_truncation",
        //    "gene_variant",
        //    "intergenic_region",
        //    "intragenic_variant",
        //    "intron_variant",
        //    "mature_miRNA_variant",
        //    "miRNA",
        //    "NMD_transcript_variant",
        //    "non_coding_transcript_exon_variant",
        //    "non_coding_transcript_variant",
        //    "regulatory_region_amplification",
        //    "regulatory_region_variant",
        //    "TF_binding_site_variant",
        //    "TFBS_amplification",
        //    "transcript_amplification",
        //    "transcript_variant",
        //    "upstream_gene_variant"
        //};

        private string[] BadTranscriptWarnings = new string[]
        {
            "WARNING_TRANSCRIPT_INCOMPLETE",
            "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS",
            "WARNING_TRANSCRIPT_NO_STOP_CODON",
            "WARNING_TRANSCRIPT_NO_START_CODON"
        };


        /// <summary>
        /// It looks like WARNING_TRANSCRIPT_INCOMPLETE, WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS,
        /// WARNING_TRANSCRIPT_NO_STOP_CODON, and WARNING_TRANSCRIPT_NO_START_CODON are relevant to this program.
        ///
        /// These are the ones that I shouldn't be translating.
        ///
        /// Could also be used for error messages regarding certain transcripts.
        /// </summary>
        public Dictionary<string, string> SnpEffWarningDescriptions = new Dictionary<string, string>
        {
            { "ERROR_CHROMOSOME_NOT_FOUND", "Chromosome does not exists in reference genome database." },
            { "ERROR_OUT_OF_CHROMOSOME_RANGE", "The variant’s genomic coordinate is greater than chromosome's length." },
            { "WARNING_REF_DOES_NOT_MATCH_GENOME", "‘REF’ in VCF does not match the reference genome." },
            { "WARNING_SEQUENCE_NOT_AVAILABLE", "Reference sequence is not available." },
            { "WARNING_TRANSCRIPT_INCOMPLETE", "Transcript length not multiple of 3 (likely incomplete in reference)." },
            { "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS", "Transcript has ≥2 internal STOP codons (possible reference error)." },
            { "WARNING_TRANSCRIPT_NO_START_CODON", "Transcript lacks START codon (possible reference error)." },
            { "WARNING_TRANSCRIPT_NO_STOP_CODON", "Transcript lacks STOP codon (possible reference error)." },
            { "INFO_REALIGN_3_PRIME", "Variant realigned to most 3′ position (HGVS compliance)." },
            { "INFO_COMPOUND_ANNOTATION", "Effect derives from compound variants." },
            { "INFO_NON_REFERENCE_ANNOTATION", "Alternative reference sequence used for annotation." },
        };
    }
}