using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;

// File location: PredictionClients/Koina/MslLibraryPredictionExtensions.cs
//
// Why an extension method rather than MslLibrary.PredictFragments():
//
//   MslLibrary  lives in Readers.
//   MslFragmentModelRouter lives in PredictionClients.
//   PredictionClients already references Readers → adding Readers → PredictionClients
//   would be circular.
//
//   Extension methods in PredictionClients can extend MslLibrary (a Readers type)
//   without modifying Readers at all. Callers simply add:
//
//       using PredictionClients.Koina;
//
//   and then write:
//
//       int updated = library.PredictFragments();

namespace PredictionClients.Koina;

/// <summary>
/// Extension methods that add Koina fragment-intensity prediction capabilities to
/// <see cref="MslLibrary"/> without requiring <c>Readers</c> to reference
/// <c>PredictionClients</c>.
/// </summary>
public static class MslLibraryPredictionExtensions
{
	/// <summary>
	/// Predicts fragment ions for all entries in <paramref name="library"/> whose
	/// <see cref="MslLibraryEntry.Source"/> is <see cref="MslFormat.SourceType.Predicted"/>
	/// and whose <see cref="MslLibraryEntry.Fragments"/> list is empty, then writes the
	/// predicted fragments back into those entries in-place.
	///
	/// <para>
	/// Routing and model selection are delegated to <see cref="MslFragmentModelRouter"/>.
	/// Each entry's stored <see cref="MslLibraryEntry.DissociationType"/> and
	/// <see cref="MslLibraryEntry.Nce"/> are used to choose the correct Koina model.
	/// Entries with <see cref="DissociationType.Unknown"/> / <c>Nce == 0</c> fall back to
	/// <see cref="MslFragmentModelRouter.DefaultDissociationType"/> /
	/// <see cref="MslFragmentModelRouter.DefaultNce"/> (HCD-28).
	/// </para>
	///
	/// <para>
	/// Entries that cannot be routed (unsupported dissociation type) are silently skipped.
	/// </para>
	///
	/// <para>
	/// <b>Modifies entries in-place.</b> Fragment ions are written directly into
	/// <see cref="MslLibraryEntry.Fragments"/> for each matched entry.
	/// </para>
	/// </summary>
	/// <param name="library">The library whose entries should be populated. Must not be null.</param>
	/// <returns>
	///   The number of entries whose <see cref="MslLibraryEntry.Fragments"/> list was
	///   successfully populated by this call.
	/// </returns>
	/// <exception cref="ObjectDisposedException">
	///   Thrown when <paramref name="library"/> has already been disposed.
	/// </exception>
	public static int PredictFragments(this MslLibrary library)
	{
		ArgumentNullException.ThrowIfNull(library);

		// Collect only Predicted-source entries with empty fragment lists
		var eligible = library
			.GetAllEntries(includeDecoys: true)
			.Where(e => e.Source == MslFormat.SourceType.Predicted && e.Fragments.Count == 0)
			.ToList();

		if (eligible.Count == 0)
			return 0;

		Dictionary<string, LibrarySpectrum> predicted =
			MslFragmentModelRouter.Predict(eligible, out _);

		if (predicted.Count == 0)
			return 0;

		var entryByKey = eligible.ToDictionary(e => e.LookupKey);
		int updatedCount = 0;

		foreach (var (key, spectrum) in predicted)
		{
			if (!entryByKey.TryGetValue(key, out MslLibraryEntry? entry))
				continue;

			// Convert LibrarySpectrum.MatchedFragmentIons → List<MslFragmentIon>
			entry.Fragments = spectrum.MatchedFragmentIons
				.Select(ion => new MslFragmentIon
				{
					Mz = (float)ion.Mz,
					Intensity = (float)ion.Intensity,
					ProductType = ion.NeutralTheoreticalProduct.ProductType,
					SecondaryProductType = null,   // Prosit does not produce internal fragments
					SecondaryFragmentNumber = 0,
					FragmentNumber = ion.NeutralTheoreticalProduct.FragmentNumber,
					ResiduePosition = ion.NeutralTheoreticalProduct.ResiduePosition,
					Charge = ion.Charge,
					NeutralLoss = ion.NeutralTheoreticalProduct.NeutralLoss,
					ExcludeFromQuant = false
				})
				.ToList();

			entry.Source = MslFormat.SourceType.Predicted;
			updatedCount++;
		}

		return updatedCount;
	}
}