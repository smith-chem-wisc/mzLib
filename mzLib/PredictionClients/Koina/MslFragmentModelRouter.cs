using MassSpectrometry;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;
using Readers.SpectralLibrary;

// File location: PredictionClients/Koina/MslFragmentModelRouter.cs
//
// Placed in PredictionClients because:
//   - Readers already imports PredictionClients types (FragmentIntensityModel, etc.)
//   - Adding PredictionClients as a reference in Readers would create a circular dependency
//   - PredictionClients already references Readers.SpectralLibrary (see FragmentIntensityModel.cs)
//     so MslLibraryEntry / MslFragmentIon are already reachable here

namespace PredictionClients.Koina;

/// <summary>
/// Routes a list of <see cref="MslLibraryEntry"/> objects to the appropriate Koina
/// fragment intensity prediction model based on each entry's stored
/// <see cref="MslLibraryEntry.DissociationType"/> and <see cref="MslLibraryEntry.Nce"/> values.
///
/// <para>
/// Grouping strategy: entries are grouped by the (DissociationType, Nce) pair.  Each unique
/// group is dispatched to one <see cref="FragmentIntensityModel"/> instance, allowing a single
/// call to predict fragments for a mixed library containing entries generated at different
/// collision energies or with different fragmentation methods.
/// </para>
///
/// <para>
/// If an entry's DissociationType is <see cref="DissociationType.Unknown"/> and Nce is 0,
/// the router substitutes <see cref="DefaultDissociationType"/> and <see cref="DefaultNce"/>.
/// Entries whose combination has no registered model are skipped; a non-null
/// <c>warnings</c> string is produced rather than throwing an exception.
/// </para>
///
/// <para>
/// <b>Extending the router:</b> add one entry to the private <c>ModelFactory</c> dictionary.
/// No other code changes are required.
/// </para>
/// </summary>
public static class MslFragmentModelRouter
{
	// ── Public constants ──────────────────────────────────────────────────────

	/// <summary>Default dissociation type used when an entry has DissociationType.Unknown.</summary>
	public const DissociationType DefaultDissociationType = DissociationType.HCD;

	/// <summary>Default NCE used when an entry has Nce == 0.</summary>
	public const int DefaultNce = 28;

	// ── Model registration table ──────────────────────────────────────────────
	//
	// Each factory produces a fresh, configured FragmentIntensityModel instance.
	// Inputs are supplied separately via model.Predict(modelInputs).
	//
	// To add a new model: one new line here, nothing else.

	private static readonly Dictionary<DissociationType, Func<FragmentIntensityModel>>
		ModelFactory = new()
		{
			[DissociationType.HCD] = () => new Prosit2020IntensityHCD(
				modHandlingMode: IncompatibleModHandlingMode.RemoveIncompatibleMods,
				parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull,
				fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence),

			// Future models — one line each:
			// [DissociationType.CID]   = () => new SomeCIDModel(...),
			// [DissociationType.EThcD] = () => new SomeEThcDModel(...),
		};

	// ── Public API ────────────────────────────────────────────────────────────

	/// <summary>
	/// Groups <paramref name="entries"/> by (DissociationType, Nce), instantiates the
	/// appropriate Koina model for each group, runs inference, and returns a dictionary
	/// mapping each entry's <see cref="MslLibraryEntry.LookupKey"/> to its predicted
	/// <see cref="LibrarySpectrum"/>.
	///
	/// <para>
	/// Entries with unsupported (DissociationType, Nce) combinations are skipped.
	/// A non-null <paramref name="warnings"/> string is produced when any entries are
	/// skipped or when the model reports internal issues (e.g. duplicate spectra).
	/// </para>
	/// </summary>
	/// <param name="entries">Library entries to predict fragments for. Must not be null.</param>
	/// <param name="warnings">
	///   Human-readable description of skipped entries or model-level warnings.
	///   <see langword="null"/> when everything succeeded without issues.
	/// </param>
	/// <returns>
	///   Dictionary from <see cref="MslLibraryEntry.LookupKey"/> to predicted
	///   <see cref="LibrarySpectrum"/>. Only successfully predicted entries are present.
	/// </returns>
	/// <exception cref="ArgumentNullException"><paramref name="entries"/> is null.</exception>
	public static Dictionary<string, LibrarySpectrum> Predict(
		IReadOnlyList<MslLibraryEntry> entries,
		out string? warnings)
	{
		ArgumentNullException.ThrowIfNull(entries);

		var results = new Dictionary<string, LibrarySpectrum>(entries.Count);
		var warnParts = new List<string>();

		if (entries.Count == 0)
		{
			warnings = null;
			return results;
		}

		// Group entries by normalised (DissociationType, Nce) key
		var groups = entries
			.GroupBy(e => NormaliseKey(e.DissociationType, e.Nce))
			.ToList();

		foreach (var group in groups)
		{
			var (dt, nce) = group.Key;
			var groupList = group.ToList();

			if (!ModelFactory.TryGetValue(dt, out var factory))
			{
				warnParts.Add(
					$"No model registered for DissociationType={dt} (NCE={nce}). " +
					$"{groupList.Count} entr{(groupList.Count == 1 ? "y" : "ies")} skipped.");
				continue;
			}

			// Build inputs — CollisionEnergy is int? in the new API
			var modelInputs = groupList
				.Select(e => new FragmentIntensityPredictionInput(
					FullSequence: e.ModifiedSequence,
					PrecursorCharge: e.Charge,
					CollisionEnergy: nce,
					InstrumentType: null,
					FragmentationType: null))
				.ToList();

			// Construct model
			FragmentIntensityModel model;
			try
			{
				model = factory();
			}
			catch (Exception ex)
			{
				warnParts.Add(
					$"Model construction failed for DissociationType={dt} NCE={nce}: {ex.Message}. " +
					$"{groupList.Count} entr{(groupList.Count == 1 ? "y" : "ies")} skipped.");
				continue;
			}

			// Run prediction
			List<PeptideFragmentIntensityPrediction> predictions;
			try
			{
				predictions = model.Predict(modelInputs);
			}
			catch (Exception ex)
			{
				warnParts.Add(
					$"Prediction failed for DissociationType={dt} NCE={nce}: {ex.Message}. " +
					$"{groupList.Count} entr{(groupList.Count == 1 ? "y" : "ies")} skipped.");
				continue;
			}

			// Collect per-entry prediction warnings
			foreach (var pred in predictions.Where(p => p.Warning != null))
				warnParts.Add($"[{pred.FullSequence}/{pred.PrecursorCharge}] {pred.Warning!.Message}");

			// Build the aligned-RT array required by GenerateLibrarySpectraFromPredictions.
			// One slot per original model input; null where iRT is unknown (0.0).
			double?[] alignedRts = BuildAlignedRts(groupList);

			// Generate LibrarySpectrum list from predictions
			List<LibrarySpectrum> spectra;
			try
			{
				spectra = model.GenerateLibrarySpectraFromPredictions(
					alignedRetentionTimes: alignedRts,
					warning: out var libWarning);

				if (libWarning is not null)
					warnParts.Add($"[DissociationType={dt} NCE={nce}] {libWarning.Message}");
			}
			catch (Exception ex)
			{
				warnParts.Add(
					$"Spectral library generation failed for DissociationType={dt} NCE={nce}: {ex.Message}.");
				continue;
			}

			// LibrarySpectrum.Name == "ModifiedSequence/ChargeState" == MslLibraryEntry.LookupKey
			foreach (var spectrum in spectra)
				results[spectrum.Name] = spectrum;
		}

		warnings = warnParts.Count > 0 ? string.Join(Environment.NewLine, warnParts) : null;
		return results;
	}

	/// <summary>
	/// Returns the concrete model <see cref="Type"/> that would be instantiated for the given
	/// dissociation type and NCE (after applying the Unknown/0 fallback rule).
	/// Returns <see langword="null"/> when no model is registered.
	/// </summary>
	public static Type? ResolveModelType(DissociationType dissociationType, int nce)
	{
		var (normDt, _) = NormaliseKey(dissociationType, nce);
		return ModelFactory.TryGetValue(normDt, out var factory)
			? factory().GetType()
			: null;
	}

	/// <summary>
	/// Returns <see langword="true"/> when the given combination has a registered model
	/// (after applying the Unknown/0 fallback rule).
	/// </summary>
	public static bool IsSupported(DissociationType dissociationType, int nce)
	{
		var (normDt, _) = NormaliseKey(dissociationType, nce);
		return ModelFactory.ContainsKey(normDt);
	}

	// ── Private helpers ───────────────────────────────────────────────────────

	private static (DissociationType dt, int nce) NormaliseKey(DissociationType dt, int nce)
	{
		bool dtUnknown = dt == DissociationType.Unknown;
		bool nceUnknown = nce == 0;

		if (dtUnknown && nceUnknown) return (DefaultDissociationType, DefaultNce);
		if (dtUnknown) return (DefaultDissociationType, nce);
		if (nceUnknown) return (dt, DefaultNce);
		return (dt, nce);
	}

	/// <summary>
	/// Builds the aligned-RT array for <see cref="FragmentIntensityModel.GenerateLibrarySpectraFromPredictions"/>.
	/// One slot per entry; null where iRT is 0.0 (unknown).
	/// </summary>
	private static double?[] BuildAlignedRts(IReadOnlyList<MslLibraryEntry> groupEntries)
	{
		var rts = new double?[groupEntries.Count];
		for (int i = 0; i < groupEntries.Count; i++)
		{
			double irt = groupEntries[i].Irt;
			rts[i] = irt == 0.0 ? null : irt;
		}
		return rts;
	}
}