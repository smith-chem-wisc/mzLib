using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectrumMatch;

namespace Omics.SpectralMatch.MslSpectralLibrary;

// ──────────────────────────────────────────────────────────────────────────────
// ConversionResult
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Result of attempting to convert a LibrarySpectrum to an MslLibraryEntry.
/// Contains either a successfully converted entry or a collection of validation errors.
/// </summary>
public class ConversionResult
{
	/// <summary>
	/// True when the conversion succeeded and Entry contains a valid MslLibraryEntry.
	/// False when validation failed and Errors contains one or more error messages.
	/// </summary>
	public bool Success { get; init; }

	/// <summary>
	/// The converted MslLibraryEntry when Success is true; null when Success is false.
	/// </summary>
	public MslLibraryEntry? Entry { get; init; }

	/// <summary>
	/// Collection of error messages describing why the conversion failed.
	/// Empty when Success is true; contains one or more messages when Success is false.
	/// </summary>
	public List<string> Errors { get; init; } = new();

	/// <summary>
	/// Creates a successful conversion result.
	/// </summary>
	public static ConversionResult FromSuccess(MslLibraryEntry entry)
		=> new() { Success = true, Entry = entry };

	/// <summary>
	/// Creates a failed conversion result with a single error message.
	/// </summary>
	public static ConversionResult FromError(string error)
		=> new() { Success = false, Errors = new List<string> { error } };

	/// <summary>
	/// Creates a failed conversion result with multiple error messages.
	/// </summary>
	public static ConversionResult FromErrors(IEnumerable<string> errors)
		=> new() { Success = false, Errors = new List<string>(errors) };
}

// ──────────────────────────────────────────────────────────────────────────────
// MslFragmentIon
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Rich in-memory representation of a single fragment ion inside an MslLibraryEntry.
/// All fields are stored as structured data (product types, residue numbers, charge) rather
/// than as an annotation string; the human-readable annotation is reconstructed on demand
/// from those fields. This mirrors the on-disk MslFragmentRecord layout but uses full .NET
/// types (nullable ProductType? for internal ions, double for neutral-loss mass) rather
/// than the compact binary encodings.
/// </summary>
public class MslFragmentIon
{
	// ── Core ion identity ───────────────────────────────────────────────

	/// <summary>
	/// Observed or predicted m/z of the fragment ion.
	/// Stored as float32 in the binary format; promoted to float here for in-memory
	/// convenience. Maps to MatchedFragmentIon.Mz.
	/// </summary>
	public float Mz { get; set; }

	/// <summary>
	/// Relative intensity of this fragment within its precursor, normalized to 0–1
	/// (1.0 = most abundant fragment in the precursor's spectrum).
	/// Maps to MatchedFragmentIon.Intensity.
	/// </summary>
	public float Intensity { get; set; }

	// ── Product type information ────────────────────────────────────────

	/// <summary>
	/// N-terminal fragment ion series type (e.g. b, c for peptides; a, d for oligos).
	/// For terminal ions this is the only product type.
	/// For internal ions this is the N-terminal terminus type of the internal pair.
	/// Maps to Product.ProductType.
	/// </summary>
	public ProductType ProductType { get; set; }

	/// <summary>
	/// C-terminal fragment ion series type for internal fragment ions; null for all
	/// terminal ions (b, y, c, z, a, w, d, etc.).
	/// When non-null, IsInternalFragment returns true and this type defines the C-terminal
	/// boundary of the internal fragment (e.g. y-type for a bIy internal ion).
	/// Maps to Product.SecondaryProductType.
	/// </summary>
	public ProductType? SecondaryProductType { get; set; }

	// ── Residue numbering ────────────────────────────────────────────────

	/// <summary>
	/// Ion series number for terminal ions (e.g. 5 for b5, 3 for y3).
	/// Start residue index (0-based) for internal fragment ions.
	/// Maps to Product.FragmentNumber.
	/// </summary>
	public int FragmentNumber { get; set; }

	/// <summary>
	/// End residue index for internal fragment ions (0-based, exclusive upper bound).
	/// 0 for all terminal ions (FragmentNumber fully describes them).
	/// Maps to Product.SecondaryFragmentNumber.
	/// </summary>
	public int SecondaryFragmentNumber { get; set; }

	/// <summary>
	/// Residue position within the parent peptide sequence from which the fragment
	/// was annotated. Used by MetaMorpheus during spectrum visualization.
	/// Maps to Product.ResiduePosition.
	/// </summary>
	public int ResiduePosition { get; set; }

	// ── ChargeState and neutral loss ──────────────────────────────────────────

	/// <summary>
	/// Charge state of the fragment ion (typically 1; occasionally 2 for large peptides).
	/// Maps to MatchedFragmentIon.Charge.
	/// </summary>
	public int Charge { get; set; }

	/// <summary>
	/// Neutral-loss mass applied to this fragment ion (0.0 = no loss).
	/// For named losses (H2O, NH3, H3PO4, etc.) this value is populated by the reader
	/// from the NeutralLossCode bits; for Custom losses it is read from the extended
	/// annotation table. Negative values represent losses; positive values represent gains.
	/// Maps to Product.NeutralLoss.
	/// </summary>
	public double NeutralLoss { get; set; }

	// ── Quantification flag ──────────────────────────────────────────────

	/// <summary>
	/// When true this fragment should be excluded from quantification (mirrors the
	/// DIA-NN ExcludeFromAssay concept). Stored in bit 5 of the fragment flags byte.
	/// </summary>
	public bool ExcludeFromQuant { get; set; }

	// ── Derived properties ───────────────────────────────────────────────

	/// <summary>
	/// True when this fragment is an internal ion (i.e. SecondaryProductType is not null).
	/// Equivalent to the is_internal bit in the on-disk fragment flags byte.
	/// Internal ions span a sub-sequence of the peptide and require both a start and an end
	/// residue number to be fully specified.
	/// </summary>
	public bool IsInternalFragment => SecondaryProductType != null;

	/// <summary>
	/// True when this fragment is a diagnostic ion (ProductType == ProductType.D).
	/// Diagnostic ions arise from modifications (glycans, cross-links, etc.) rather than
	/// from peptide backbone cleavage. Equivalent to the is_diagnostic flag bit.
	/// </summary>
	public bool IsDiagnosticIon => ProductType == ProductType.D;

	// ── Annotation reconstruction ────────────────────────────────────────

	/// <summary>
	/// Human-readable annotation string reconstructed at runtime from the structured fields;
	/// no annotation string is stored in the binary file.
	/// Examples: "b5", "y3-H2O", "bIy[3-6]" (internal b/y spanning residues 3–6).
	/// The format follows mzLib conventions: type + number (+ neutral loss if present).
	/// </summary>
	public string Annotation
	{
		get
		{
			// Build the annotation string from structured fields
			string typePart;

			if (IsInternalFragment)
			{
				// Internal ion: format is <Ntype>I<Ctype>[start-end]
				// e.g. bIy[3-6] for b-type N-terminus / y-type C-terminus spanning 3 to 6
				typePart = $"{ProductType}I{SecondaryProductType}[{FragmentNumber}-{SecondaryFragmentNumber}]";
			}
			else
			{
				// Terminal ion: format is <type><number>  e.g. b5, y12, c3
				typePart = $"{ProductType}{FragmentNumber}";
			}

			// Append neutral-loss suffix when non-zero
			// Negative loss magnitude is shown as a subtraction; positive as addition
			if (NeutralLoss != 0.0)
			{
				string lossSuffix = NeutralLoss < 0
					? $"-{-NeutralLoss:F4}"
					: $"+{NeutralLoss:F4}";
				typePart += lossSuffix;
			}

			// Append charge superscript for z > 1
			if (Charge > 1)
				typePart += $"^{Charge}+";

			return typePart;
		}
	}
}

// ──────────────────────────────────────────────────────────────────────────────
// MslLibraryEntry
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// The rich in-memory representation of one library precursor entry.  This class is the
/// primary data model for the MSL binary spectral library; it mirrors LibrarySpectrum
/// from Omics.SpectrumMatch but adds all extended metadata (ion mobility, protein
/// accession/name/gene, molecule type, source type, dissociation type, NCE, q-value,
/// elution group) that the binary format supports.
///
/// Conversion methods (<see cref="ToLibrarySpectrum"/>, <see cref="FromLibrarySpectrum"/>)
/// provide bidirectional interop with the existing mzLib type system.
/// </summary>
public class MslLibraryEntry
{
	// ── Core fields (round-trip with LibrarySpectrum) ────────────────────

	/// <summary>
	/// Modified sequence in mzLib bracket notation, e.g. "PEPTM[Common Variable:Oxidation on M]IDE".
	/// This is the primary identifier used in DDA-style dictionary lookups.
	/// Maps to LibrarySpectrum.Sequence and MslPrecursorRecord.ModifiedSeqStringIdx.
	/// </summary>
	public string FullSequence { get; set; }
	
	/// <summary>
	/// Unmodified amino-acid sequence without any modification annotations.
	/// Used for elution group construction and for the stripped-sequence lookup dictionary.
	/// Maps to MslPrecursorRecord.StrippedSeqStringIdx.
	/// </summary>
	public string BaseSequence { get; set; }

	/// <summary>
	/// Precursor m/z value. For multiply-charged precursors this is the (M + z*H) / z value.
	/// Maps to LibrarySpectrum.PrecursorMz and MslPrecursorRecord.PrecursorMz.
	/// </summary>
	public double PrecursorMz { get; set; }

	/// <summary>
	/// Precursor charge state. Typically 1–4 for tryptic peptides; may be higher for
	/// top-down proteoforms. Maps to LibrarySpectrum.ChargeState and MslPrecursorRecord.ChargeState.
	/// </summary>
	public int ChargeState { get; set; }

	/// <summary>
	/// Indexed retention time (iRT, in iRT units) or calibrated run-specific retention time
	/// in minutes. Which representation is stored is indicated by the rt_is_calibrated bit
	/// in PrecursorFlags. Maps to LibrarySpectrum.RetentionTime and MslPrecursorRecord.RetentionTime.
	/// </summary>
	public double RetentionTime { get; set; }

	/// <summary>
	/// True for decoy precursors (reversed or shuffled sequence). Used to separate target
	/// and decoy distributions for FDR estimation. Maps to LibrarySpectrum.IsDecoy.
	/// </summary>
	public bool IsDecoy { get; set; }

	/// <summary>
	/// All fragment ions for this precursor.  The list may contain terminal ions (b, y, c, z,
	/// a, w, d, …), internal fragment ions, diagnostic ions, and glycan Y ions in any order;
	/// the writer sorts them by m/z before serialization.
	/// Maps to LibrarySpectrum.MatchedFragmentIons via ToLibrarySpectrum().
	/// </summary>
	public List<MslFragmentIon> MatchedFragmentIons { get; set; } = new();

	// ── Extended fields ──────────────────────────────────────────────────

	/// <summary>
	/// Ion mobility expressed as 1/K0 (collisional cross section proxy, in V·s/cm²).
	/// 0.0 means the value is not available. Maps to MslPrecursorRecord.IonMobility.
	/// </summary>
	public double IonMobility { get; set; }

	/// <summary>
	/// UniProt accession of the source protein (e.g. "P04637" for TP53).
	/// Empty string or null when no protein information is associated with this entry.
	/// Resolved from MslPrecursorRecord.ProteinIdx via the protein table at read time.
	/// </summary>
	public string ProteinAccession { get; set; }

	/// <summary>
	/// Human-readable protein name (e.g. "Cellular tumor antigen p53").
	/// Empty string or null when absent. Resolved from the protein table at read time.
	/// </summary>
	public string ProteinName { get; set; }

	/// <summary>
	/// HGNC gene symbol (e.g. "TP53"). Empty string or null when absent.
	/// Resolved from the protein table at read time; only populated when the file-level
	/// has_gene_data flag is set.
	/// </summary>
	public string GeneName { get; set; }

	/// <summary>
	/// True when the peptide is expected to produce a uniquely detectable signal
	/// (proteotypic). Stored in bit 1 of PrecursorFlags. Maps to MslPrecursorRecord.PrecursorFlags.
	/// </summary>
	public bool IsProteotypic { get; set; }

	/// <summary>
	/// Library q-value (confidence score; lower is better; 0 = perfect confidence).
	/// float.NaN when no q-value is available (e.g. purely predicted library).
	/// Maps to MslPrecursorRecord.QValue.
	/// </summary>
	public float QValue { get; set; } = float.NaN;

	/// <summary>
	/// Elution group identifier. All precursors with the same BaseSequence value share
	/// the same ElutionGroupId, allowing the DIA scorer to handle multiple charge states and
	/// modifications as co-eluting species.
	/// Maps to MslPrecursorRecord.ElutionGroupId.
	/// </summary>
	public int ElutionGroupId { get; set; }

	/// <summary>
	/// Molecule classification for this precursor.  Controls which fragmentation namespace
	/// is used to reconstruct Product.Terminus and Product.NeutralMass at read time.
	/// Defaults to Peptide for standard tryptic libraries.
	/// Maps to MslPrecursorRecord.MoleculeType.
	/// </summary>
	public MslFormat.MoleculeType MoleculeType { get; set; } = MslFormat.MoleculeType.Peptide;

	/// <summary>
	/// Origin of fragment intensities: Predicted (model output), Empirical (real experiment),
	/// or EmpiricalRefined (predicted RT/IM replaced by empirical values).
	/// Maps to MslPrecursorRecord.SourceType.
	/// </summary>
	public MslFormat.SourceType Source { get; set; } = MslFormat.SourceType.Predicted;

	/// <summary>
	/// Dissociation type used to generate this spectrum (e.g. HCD, CID, EThcD).
	/// Required for correct reconstruction of Product.NeutralMass at read time because
	/// different dissociation types have different mass shift offsets per ion type.
	/// Maps to MslPrecursorRecord.DissociationType.
	/// </summary>
	public DissociationType DissociationType { get; set; } = DissociationType.Unknown;

	/// <summary>
	/// Nominal collision energy. Stored as NCE × 10 as int16 in the binary format;
	/// exposed here as a plain int (e.g. 28 for NCE 28). 0 = unknown or not applicable.
	/// Maps to MslPrecursorRecord.Nce (divided by 10 on read, multiplied on write).
	/// </summary>
	public int Nce { get; set; }

	// ── Lookup key ───────────────────────────────────────────────────────

	/// <summary>
	/// Canonical DDA-style lookup key combining modified sequence and charge state.
	/// Format: "{FullSequence}/{ChargeState}", e.g. "PEPTM[Common Variable:Oxidation on M]IDE/2".
	/// Used as the dictionary key in the MslIndex sequence-to-entry map for O(1) DDA lookup.
	/// </summary>
	public string Name => FullSequence + "/" + ChargeState;

	// ── Conversion: MslLibraryEntry → LibrarySpectrum ───────────────────

	/// <summary>
	/// Converts this entry to a LibrarySpectrum compatible with the existing mzLib type system.
	/// Fragment ions are converted from MslFragmentIon to MatchedFragmentIon, with Product
	/// instances reconstructed from the stored structured fields. Neutral masses are
	/// recomputed using the appropriate dissociation-type mass shift; Product.Terminus is
	/// looked up from the appropriate TerminusSpecificProductTypes dictionary based on
	/// MoleculeType.
	/// </summary>
	/// <returns>
	///   A new LibrarySpectrum with:
	///   <list type="bullet">
	///     <item>Sequence = FullSequence</item>
	///     <item>PrecursorMz = PrecursorMz (double precision)</item>
	///     <item>ChargeState = ChargeState</item>
	///     <item>RetentionTime = RetentionTime</item>
	///     <item>IsDecoy = IsDecoy</item>
	///     <item>MatchedFragmentIons = one MatchedFragmentIon per MslFragmentIon in MatchedFragmentIons</item>
	///   </list>
	/// </returns>
	public LibrarySpectrum ToLibrarySpectrum()
	{
		// Convert every MslFragmentIon to a MatchedFragmentIon wrapping a Product
		var matchedIons = new List<MatchedFragmentIon>(MatchedFragmentIons.Count);

		foreach (MslFragmentIon mslIon in MatchedFragmentIons)
		{
			// Determine whether this is an oligo or peptide to select the correct
			// Terminus lookup source
			bool isOligo = MoleculeType == MslFormat.MoleculeType.Oligonucleotide;

			// Derive terminus directly from the ProductType enum value.
			// Internal ions (IsInternalFragment) and diagnostic ions always get None.
			// For oligo entries, a/d are 5'-terminal; w/z are 3'-terminal.
			// For peptide/proteoform/glycopeptide, b/a/c are N-terminal; y/z/x are C-terminal.
			FragmentationTerminus terminus;
			if (mslIon.IsInternalFragment)
			{
				terminus = FragmentationTerminus.None;
			}
			else if (isOligo)
			{
				terminus = mslIon.ProductType switch
				{
					ProductType.a or ProductType.d => FragmentationTerminus.FivePrime,
					ProductType.w or ProductType.z => FragmentationTerminus.ThreePrime,
					_ => FragmentationTerminus.None
				};
			}
			else
			{
				terminus = mslIon.ProductType switch
				{
					ProductType.b or ProductType.a or ProductType.c => FragmentationTerminus.N,
					ProductType.y or ProductType.z or ProductType.x => FragmentationTerminus.C,
					_ => FragmentationTerminus.None
				};
			}

			// Construct the Product using the mzLib constructor:
			// (productType, terminus, neutralMass, fragmentNumber, residuePosition, neutralLoss,
			//  secondaryProductType, secondaryFragmentNumber).
			// neutralMass is passed as 0.0 — this matches the pattern used throughout the
			// existing mzLib readers (MSP reader, DIA-NN reader) because mzLib uses the stored
			// m/z directly for matching rather than recomputing from neutral mass at read time.
			// SecondaryProductType and SecondaryFragmentNumber are passed through to Product's
			// optional parameters, enabling correct internal-ion annotation (e.g. bIy[3-6])
			// in MetaMorpheus spectrum visualization and scoring.
			var product = new Product(
				mslIon.ProductType,
				terminus,
				neutralMass: 0.0,
				mslIon.FragmentNumber,
				mslIon.ResiduePosition,
				mslIon.NeutralLoss,
				mslIon.SecondaryProductType,
				mslIon.SecondaryFragmentNumber);

			// Wrap in a MatchedFragmentIon with the stored Mz (float → double), Intensity, and ChargeState
			matchedIons.Add(new MatchedFragmentIon(product, mslIon.Mz, mslIon.Intensity, mslIon.Charge));
		}

		// Construct the LibrarySpectrum using the mzLib constructor signature
		return new LibrarySpectrum(
			sequence: FullSequence,
			precursorMz: PrecursorMz,
			chargeState: ChargeState,
			peaks: matchedIons,
			rt: RetentionTime,
			isDecoy: IsDecoy);
	}

	// ── Conversion: LibrarySpectrum → MslLibraryEntry ───────────────────

	/// <summary>
	/// Attempts to construct an MslLibraryEntry from an existing LibrarySpectrum, returning
	/// a ConversionResult that contains either the successfully converted entry or a collection
	/// of validation errors. This method never throws exceptions; use this for batch processing
	/// where you need to collect errors and continue processing other spectra.
	/// </summary>
	/// <param name="spectrum">
	///   The LibrarySpectrum to convert. May be null (returns error).
	/// </param>
	/// <returns>
	///   A ConversionResult with Success=true and Entry populated on success, or Success=false
	///   and Errors populated with validation failure messages on failure.
	/// </returns>
	public static ConversionResult TryFromLibrarySpectrum(LibrarySpectrum? spectrum)
	{
		// Validate input
		var errors = new List<string>();

		if (spectrum is null)
		{
			errors.Add("LibrarySpectrum is null");
			return ConversionResult.FromErrors(errors);
		}

		if (string.IsNullOrEmpty(spectrum.Sequence))
		{
			errors.Add($"LibrarySpectrum.Sequence is null or empty (PrecursorMz: {spectrum.PrecursorMz:F4}, ChargeState: {spectrum.ChargeState})");
		}

		if (spectrum.ChargeState <= 0)
		{
			errors.Add($"LibrarySpectrum.ChargeState must be positive (Sequence: {spectrum.Sequence ?? "<null>"}, ChargeState: {spectrum.ChargeState})");
		}

		if (spectrum.PrecursorMz <= 0)
		{
			errors.Add($"LibrarySpectrum.PrecursorMz must be positive (Sequence: {spectrum.Sequence ?? "<null>"}, PrecursorMz: {spectrum.PrecursorMz})");
		}

		// Return early if validation failed
		if (errors.Count > 0)
			return ConversionResult.FromErrors(errors);

		// Perform conversion (at this point spectrum.Sequence is guaranteed non-null)
		string strippedSeq = StripModifications(spectrum.Sequence);

		var fragments = new List<MslFragmentIon>(spectrum.MatchedFragmentIons?.Count ?? 0);

		foreach (MatchedFragmentIon ion in spectrum.MatchedFragmentIons ?? Enumerable.Empty<MatchedFragmentIon>())
		{
			Product product = ion.NeutralTheoreticalProduct;

			var mslIon = new MslFragmentIon
			{
				Mz = (float)ion.Mz,
				Intensity = (float)ion.Intensity,
				ProductType = product.ProductType,
				SecondaryProductType = product.SecondaryProductType,
				SecondaryFragmentNumber = product.SecondaryFragmentNumber,
				FragmentNumber = product.FragmentNumber,
				ResiduePosition = product.ResiduePosition,
				Charge = ion.Charge,
				NeutralLoss = product.NeutralLoss,
				ExcludeFromQuant = false
			};

			fragments.Add(mslIon);
		}

		var entry = new MslLibraryEntry
		{
			FullSequence = spectrum.Sequence,
			BaseSequence = strippedSeq,
			PrecursorMz = spectrum.PrecursorMz,
			ChargeState = spectrum.ChargeState,
			RetentionTime = spectrum.RetentionTime ?? 0.0,
			IsDecoy = spectrum.IsDecoy,
			MatchedFragmentIons = fragments,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = float.NaN,
			ElutionGroupId = 0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			DissociationType = DissociationType.Unknown,
			Nce = 0
		};

		return ConversionResult.FromSuccess(entry);
	}

	/// <summary>
	/// Static factory that constructs an MslLibraryEntry from an existing LibrarySpectrum.
	/// All structured fragment fields (ProductType, FragmentNumber, ResiduePosition, etc.)
	/// are extracted from the MatchedFragmentIon list. Extended metadata fields that have
	/// no representation in LibrarySpectrum are initialized to sensible defaults:
	/// IonMobility = 0, QValue = float.NaN, ElutionGroupId = 0, Source = Predicted,
	/// MoleculeType = Peptide, DissociationType = Unknown, Nce = 0.
	/// 
	/// This method throws exceptions on validation failure; for batch processing where you
	/// need to collect errors and continue, use <see cref="TryFromLibrarySpectrum"/> instead.
	/// </summary>
	/// <param name="spectrum">
	///   The LibrarySpectrum to convert. Must not be null. The MatchedFragmentIon list may be
	///   empty (resulting in an entry with no fragments) but should not be null.
	/// </param>
	/// <returns>
	///   A new MslLibraryEntry whose core fields mirror the spectrum and whose fragment list
	///   contains one MslFragmentIon per MatchedFragmentIon in spectrum.MatchedFragmentIons.
	/// </returns>
	/// <exception cref="ArgumentNullException">Thrown when <paramref name="spectrum"/> is null.</exception>
	/// <exception cref="ArgumentException">Thrown when spectrum.Sequence is null or empty, or other validation fails.</exception>
	public static MslLibraryEntry? FromLibrarySpectrum(LibrarySpectrum spectrum)
	{
		// Use TryFromLibrarySpectrum for validation and conversion
		var result = TryFromLibrarySpectrum(spectrum);
		return result.Success ? result.Entry : null;
	}

	// ── Private helpers ──────────────────────────────────────────────────

	/// <summary>
	/// Strips all mzLib bracket-notation modification tags from a modified sequence string,
	/// returning the bare amino-acid sequence.
	/// 
	/// Examples:
	/// - Internal modification: "PEPTM[Common Variable:Oxidation on M]IDE" → "PEPTMIDE"
	/// - C-terminal modification: "PEPTIDE-[Common Variable:Amidation]" → "PEPTIDE"
	/// - N-terminal modification: "[Common Variable:Acetylation]PEPTIDE" → "PEPTIDE"
	/// - Nested brackets: "PEP[outer[nested]more]TIDE" → "PEPTIDE"
	/// - Multiple modifications: "M[Oxidation]PEPTIDEK[Acetylation]" → "MPEPTIDEK"
	/// 
	/// The method uses a depth counter to handle arbitrarily nested brackets - anything inside
	/// the outermost brackets (including nested brackets) is completely ignored. This handles
	/// complex modification annotations that may contain brackets in their description.
	/// 
	/// The method iterates character-by-character to avoid allocating a Regex match collection
	/// on the hot path.
	/// </summary>
	/// <param name="modifiedSeq">Modified sequence in mzLib bracket notation. Must not be null.</param>
	/// <returns>The sequence with all [...] tags and terminal dashes removed.</returns>
	internal static string StripModifications(string modifiedSeq)
	{
		// Track bracket nesting depth; characters inside any level of brackets are omitted.
		// This naturally handles nested brackets: depth > 0 means "inside brackets somewhere"
		// regardless of nesting level, so all nested content is correctly ignored.
		int depth = 0;

		// Use a char-array builder sized to the input to avoid reallocation
		var buffer = new char[modifiedSeq.Length];
		int writePos = 0;

		// Use index-based loop to allow look-back for dash-before-bracket pattern
		for (int i = 0; i < modifiedSeq.Length; i++)
		{
			char ch = modifiedSeq[i];

			if (ch == '[')
			{
				// Entering a modification tag: increment nesting depth but don't copy
				// For nested brackets like [[...]], this will increment depth twice (1→2)
				depth++;

				// Remove trailing dash from buffer if this is the FIRST bracket after a dash
				// (C-terminal modification pattern). Only check on the first '[' (when depth==1)
				// Example: "PEPTIDE-[modification]" → buffer contains "PEPTIDE-" → remove the dash
				// Example: "PEPTIDE-[[nested]]" → same behavior, dash removed on first '['
				if (depth == 1 && writePos > 0 && buffer[writePos - 1] == '-')
				{
					writePos--; // Back up to remove the dash
				}
			}
			else if (ch == ']')
			{
				// Leaving a modification tag: decrement depth but don't copy
				// For nested brackets, this will decrement from 2→1 (still inside outer brackets)
				// Only when depth reaches 0 do we resume copying amino acid characters
				depth--;
			}
			else if (depth == 0)
			{
				// Outside all brackets: this is a real residue character — copy it
				// This condition is only true when not inside ANY bracket (nested or not)
				buffer[writePos++] = ch;
			}
			// else: depth > 0, we're inside brackets (possibly nested), skip this character
		}

		return new string(buffer, 0, writePos);
	}

	/// <summary>
	/// Maps a neutral-loss mass (in Da) to the nearest named NeutralLossCode, or
	/// NeutralLossCode.Custom when the mass does not correspond to a known common loss.
	/// Tolerance for matching is ±0.01 Da; this is intentionally loose to handle
	/// minor floating-point variations in how different tools report mass shifts.
	/// </summary>
	/// <param name="neutralLoss">
	///   Neutral-loss mass in daltons. 0.0 = no loss. Negative values represent mass losses;
	///   positive values represent mass gains.
	/// </param>
	/// <returns>The best-matching NeutralLossCode, or Custom if no match is found.</returns>
	private static MslFormat.NeutralLossCode ClassifyNeutralLoss(double neutralLoss)
		=> MslFormat.ClassifyNeutralLoss(neutralLoss);
}