using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Represents a single fragment ion in DIA-NN's data model.
    /// This is a lightweight, format-neutral representation used for TSV/Parquet/binary I/O.
    /// Converts to mzLib's MatchedFragmentIon for integration with MetaMorpheus scoring.
    /// </summary>
    public class DiaNNFragmentIon
    {
        /// <summary>Fragment m/z value</summary>
        public double Mz { get; set; }

        /// <summary>Relative intensity, typically normalized to [0, 1] with most abundant = 1.0</summary>
        public double Intensity { get; set; }

        /// <summary>Ion series type: 'b' (N-terminal) or 'y' (C-terminal)</summary>
        public char IonType { get; set; }

        /// <summary>Fragment series number (1-based, e.g., b3 → 3, y7 → 7)</summary>
        public int SeriesNumber { get; set; }

        /// <summary>Fragment ion charge state</summary>
        public int Charge { get; set; }

        /// <summary>
        /// Neutral loss type string: "noloss", "H2O", "NH3", "H3PO4", "HPO3".
        /// Matches DIA-NN's FragmentLossType TSV column.
        /// </summary>
        public string LossType { get; set; } = "noloss";

        /// <summary>
        /// Whether this fragment should be excluded from quantification assays.
        /// Corresponds to DIA-NN's ExcludeFromAssay column.
        /// </summary>
        public bool ExcludeFromAssay { get; set; } = false;
    }

    /// <summary>
    /// Rich wrapper class for a DIA-NN spectral library entry.
    /// 
    /// Contains all DIA-NN-specific metadata (ion mobility, protein groups, gene names,
    /// proteotypicity, q-values) that are not present in mzLib's LibrarySpectrum.
    /// 
    /// Design decision (from architecture doc, Section 7.3):
    /// We use a wrapper (Option B) rather than extending LibrarySpectrum to keep the
    /// core mzLib model lean. DIA-NN metadata is important for DIA workflows but
    /// irrelevant for DDA spectral matching.
    /// 
    /// Conversion:
    /// - ToLibrarySpectrum() → for MetaMorpheus scoring, spectral angle computation
    /// - FromLibrarySpectrum() → for writing mzLib data to DIA-NN formats
    /// </summary>
    public class DiaNNLibraryEntry
    {
        #region Core Spectrum Data (maps to LibrarySpectrum)

        /// <summary>
        /// Modified peptide sequence in DIA-NN format: _PEPTM[UniMod:35]IDE_
        /// Flanking underscores included. UniMod notation for modifications.
        /// </summary>
        public string ModifiedSequence { get; set; }

        /// <summary>
        /// Unmodified peptide sequence (no modifications, no underscores): PEPTMIDE
        /// </summary>
        public string StrippedSequence { get; set; }

        /// <summary>Precursor m/z value</summary>
        public double PrecursorMz { get; set; }

        /// <summary>Precursor charge state</summary>
        public int PrecursorCharge { get; set; }

        /// <summary>
        /// Retention time value. May be indexed RT (iRT) or calibrated RT in minutes,
        /// depending on the library source. Corresponds to DIA-NN's Tr_recalibrated column.
        /// </summary>
        public double RetentionTime { get; set; }

        /// <summary>Whether this entry is a decoy (reversed/shuffled sequence)</summary>
        public bool IsDecoy { get; set; }

        /// <summary>Fragment ions for this precursor</summary>
        public List<DiaNNFragmentIon> Fragments { get; set; } = new();

        #endregion

        #region DIA-NN-Specific Metadata

        /// <summary>
        /// Ion mobility value (1/K0 for timsTOF TIMS-DIA).
        /// 0.0 or NaN if not available (non-ion-mobility instruments).
        /// </summary>
        public double IonMobility { get; set; }

        /// <summary>
        /// Protein accession(s), semicolon-separated for shared peptides.
        /// Corresponds to DIA-NN's ProteinId column.
        /// </summary>
        public string ProteinId { get; set; } = string.Empty;

        /// <summary>
        /// Protein description/name(s).
        /// Corresponds to DIA-NN's ProteinName column.
        /// </summary>
        public string ProteinName { get; set; } = string.Empty;

        /// <summary>
        /// Gene name(s), semicolon-separated.
        /// Corresponds to DIA-NN's Genes column.
        /// </summary>
        public string GeneName { get; set; } = string.Empty;

        /// <summary>
        /// Whether this peptide is proteotypic (unique to one protein group).
        /// True = proteotypic, False = shared across multiple proteins.
        /// </summary>
        public bool IsProteotypic { get; set; }

        /// <summary>
        /// Library-level q-value for FDR control. Null if not available.
        /// </summary>
        public double? QValue { get; set; }

        #endregion

        #region Lookup Key

        /// <summary>
        /// Unique identifier for this entry, matching LibrarySpectrum.Name format.
        /// Uses mzLib sequence format (no underscores, descriptive mod names) + "/" + charge.
        /// </summary>
        public string Name
        {
            get
            {
                string mzLibSeq = DiaNNModificationMapping.DiaNNToMzLib(ModifiedSequence);
                return mzLibSeq + "/" + PrecursorCharge;
            }
        }

        #endregion

        #region Conversion: DiaNNLibraryEntry → LibrarySpectrum

        /// <summary>
        /// Converts this DIA-NN entry to an mzLib LibrarySpectrum for use in
        /// MetaMorpheus scoring, spectral angle computation, and other mzLib workflows.
        /// 
        /// The conversion:
        /// 1. Converts the modified sequence from DIA-NN UniMod format to mzLib descriptive format
        /// 2. Creates a Product for each fragment (ion type, terminus, fragment number, neutral loss)
        /// 3. Creates a MatchedFragmentIon for each fragment (m/z, intensity, charge)
        /// 4. Constructs a LibrarySpectrum with the converted data
        /// 
        /// DIA-NN-specific metadata (ion mobility, protein ID, gene name, etc.) is NOT carried
        /// over to LibrarySpectrum — use the wrapper directly when that metadata is needed.
        /// </summary>
        /// <returns>An mzLib LibrarySpectrum with the same spectral data</returns>
        public LibrarySpectrum ToLibrarySpectrum()
        {
            string mzLibSequence = DiaNNModificationMapping.DiaNNToMzLib(ModifiedSequence);
            string strippedSeq = StrippedSequence ?? DiaNNModificationMapping.GetStrippedSequence(ModifiedSequence);
            int peptideLength = strippedSeq.Length;

            var matchedIons = new List<MatchedFragmentIon>(Fragments.Count);

            foreach (var fragment in Fragments)
            {
                // Determine fragmentation terminus from ion type
                FragmentationTerminus terminus = fragment.IonType switch
                {
                    'b' or 'B' => FragmentationTerminus.N,
                    'y' or 'Y' => FragmentationTerminus.C,
                    _ => FragmentationTerminus.None
                };

                // Parse the ion type to mzLib ProductType
                ProductType productType = (ProductType)Enum.Parse(typeof(ProductType),
                    fragment.IonType.ToString(), ignoreCase: true);

                // Compute neutral loss mass from the DIA-NN loss type string
                double neutralLoss = 0;
                if (!string.IsNullOrEmpty(fragment.LossType) &&
                    DiaNNModificationMapping.NeutralLossNameToMass.TryGetValue(fragment.LossType, out var lossMass))
                {
                    neutralLoss = lossMass;
                }

                // Compute residue position: for b-ions it equals fragment number,
                // for y-ions it equals peptideLength - fragment number
                int residuePosition = terminus == FragmentationTerminus.N
                    ? fragment.SeriesNumber
                    : peptideLength - fragment.SeriesNumber;

                // Create the theoretical product (ion annotation)
                var product = new Product(
                    productType: productType,
                    terminus: terminus,
                    neutralMass: 0.0, // Not computed here — mzLib uses m/z directly for matching
                    fragmentNumber: fragment.SeriesNumber,
                    residuePosition: residuePosition,
                    neutralLoss: neutralLoss
                );

                // Create the matched fragment ion (observed data)
                var matchedIon = new MatchedFragmentIon(
                    neutralTheoreticalProduct: product,
                    experMz: fragment.Mz,
                    experIntensity: fragment.Intensity,
                    charge: fragment.Charge
                );

                matchedIons.Add(matchedIon);
            }

            return new LibrarySpectrum(
                sequence: mzLibSequence,
                precursorMz: PrecursorMz,
                chargeState: PrecursorCharge,
                peaks: matchedIons,
                rt: RetentionTime,
                isDecoy: IsDecoy
            );
        }

        #endregion

        #region Conversion: LibrarySpectrum → DiaNNLibraryEntry

        /// <summary>
        /// Creates a DiaNNLibraryEntry from an mzLib LibrarySpectrum.
        /// Used when writing mzLib spectral data to DIA-NN formats (TSV, Parquet, binary).
        /// 
        /// DIA-NN-specific metadata (ion mobility, protein ID, gene name, etc.) will be
        /// set to default values and can be populated separately.
        /// </summary>
        /// <param name="spectrum">Source LibrarySpectrum from mzLib</param>
        /// <returns>A DiaNNLibraryEntry with spectral data converted from mzLib format</returns>
        public static DiaNNLibraryEntry FromLibrarySpectrum(LibrarySpectrum spectrum)
        {
            string diannSequence = DiaNNModificationMapping.MzLibToDiaNN(spectrum.Sequence);
            string strippedSequence = DiaNNModificationMapping.GetStrippedSequence(spectrum.Sequence);

            var fragments = new List<DiaNNFragmentIon>(spectrum.MatchedFragmentIons.Count);

            foreach (var ion in spectrum.MatchedFragmentIons)
            {
                string ionTypeStr = ion.NeutralTheoreticalProduct.ProductType.ToString().ToLowerInvariant();
                char ionTypeChar = ionTypeStr.Length > 0 ? ionTypeStr[0] : 'y';

                string lossType = DiaNNModificationMapping.MassToNeutralLossName(
                    ion.NeutralTheoreticalProduct.NeutralLoss);

                fragments.Add(new DiaNNFragmentIon
                {
                    Mz = ion.Mz,
                    Intensity = ion.Intensity,
                    IonType = ionTypeChar,
                    SeriesNumber = ion.NeutralTheoreticalProduct.FragmentNumber,
                    Charge = ion.Charge,
                    LossType = lossType
                });
            }

            return new DiaNNLibraryEntry
            {
                ModifiedSequence = diannSequence,
                StrippedSequence = strippedSequence,
                PrecursorMz = spectrum.PrecursorMz,
                PrecursorCharge = spectrum.ChargeState,
                RetentionTime = spectrum.RetentionTime ?? 0.0,
                IsDecoy = spectrum.IsDecoy,
                Fragments = fragments,
                // DIA-NN-specific metadata defaults
                IonMobility = 0.0,
                ProteinId = string.Empty,
                ProteinName = string.Empty,
                GeneName = string.Empty,
                IsProteotypic = true,
                QValue = null,
            };
        }

        #endregion
    }
}
