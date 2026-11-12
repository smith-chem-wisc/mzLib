using Omics.Fragmentation;
using System;

namespace Omics.SpectralMatch.MslSpectralLibrary;

/// <summary>
/// One fragment ion of a <c>.msl</c> library entry, as a value.
///
/// <para>
/// It has the same fields as <see cref="MslFragmentIon"/>, decoded the same way, in a readonly
/// struct. <c>MslLibrary.ReadFragments</c> copies an entry's fragments into a caller-owned buffer
/// of these, so reading fragments allocates nothing.
/// </para>
///
/// <para>
/// Unlike the on-disk <see cref="MslFragmentRecord"/>, every field is decoded:
/// <see cref="NeutralLoss"/> is a mass, and <see cref="ResiduePosition"/> is a residue position
/// (0 for a fragment with a custom neutral loss, as in <see cref="MslFragmentIon"/>).
/// </para>
///
/// <para>
/// There is no <c>Annotation</c>: <see cref="MslFragmentIon.Annotation"/> builds that string from
/// these fields, so build an <see cref="MslFragmentIon"/> (or call <c>MslLibrary.GetEntry</c>) when
/// the text is needed.
/// </para>
/// </summary>
public readonly struct MslFragmentValue
{
	/// <summary>
	/// Creates a fragment value. Each argument sets the property of the same name.
	/// </summary>
	public MslFragmentValue(float mz, float intensity, ProductType productType,
		ProductType? secondaryProductType, int fragmentNumber, int secondaryFragmentNumber,
		int residuePosition, int charge, double neutralLoss, bool excludeFromQuant)
	{
		Mz = mz;
		Intensity = intensity;
		ProductType = productType;
		SecondaryProductType = secondaryProductType;
		FragmentNumber = fragmentNumber;
		SecondaryFragmentNumber = secondaryFragmentNumber;
		ResiduePosition = residuePosition;
		Charge = charge;
		NeutralLoss = neutralLoss;
		ExcludeFromQuant = excludeFromQuant;
	}

	/// <summary>
	/// Creates a value with the fields of <paramref name="ion"/>.
	/// </summary>
	/// <exception cref="ArgumentNullException"><paramref name="ion"/> is null.</exception>
	public MslFragmentValue(MslFragmentIon ion)
	{
		ArgumentNullException.ThrowIfNull(ion);
		Mz = ion.Mz;
		Intensity = ion.Intensity;
		ProductType = ion.ProductType;
		SecondaryProductType = ion.SecondaryProductType;
		FragmentNumber = ion.FragmentNumber;
		SecondaryFragmentNumber = ion.SecondaryFragmentNumber;
		ResiduePosition = ion.ResiduePosition;
		Charge = ion.Charge;
		NeutralLoss = ion.NeutralLoss;
		ExcludeFromQuant = ion.ExcludeFromQuant;
	}

	/// <summary>Observed or predicted m/z of the fragment ion. See <see cref="MslFragmentIon.Mz"/>.</summary>
	public float Mz { get; }

	/// <summary>
	/// Relative intensity within the precursor, 0-1 (1.0 = most abundant fragment).
	/// See <see cref="MslFragmentIon.Intensity"/>.
	/// </summary>
	public float Intensity { get; }

	/// <summary>
	/// Ion series type; for internal ions, the N-terminal type of the pair.
	/// See <see cref="MslFragmentIon.ProductType"/>.
	/// </summary>
	public ProductType ProductType { get; }

	/// <summary>
	/// C-terminal type for internal ions; null for terminal ions.
	/// See <see cref="MslFragmentIon.SecondaryProductType"/>.
	/// </summary>
	public ProductType? SecondaryProductType { get; }

	/// <summary>
	/// Ion series number for terminal ions; start residue (0-based) for internal ions.
	/// See <see cref="MslFragmentIon.FragmentNumber"/>.
	/// </summary>
	public int FragmentNumber { get; }

	/// <summary>
	/// End residue (0-based, exclusive) for internal ions; 0 for terminal ions.
	/// See <see cref="MslFragmentIon.SecondaryFragmentNumber"/>.
	/// </summary>
	public int SecondaryFragmentNumber { get; }

	/// <summary>
	/// Residue position the fragment was annotated from; 0 for a fragment with a custom neutral
	/// loss. See <see cref="MslFragmentIon.ResiduePosition"/>.
	/// </summary>
	public int ResiduePosition { get; }

	/// <summary>Charge state of the fragment ion. See <see cref="MslFragmentIon.Charge"/>.</summary>
	public int Charge { get; }

	/// <summary>
	/// Neutral-loss mass in daltons: positive for a loss, as in <c>Product.NeutralLoss</c>; 0.0 for
	/// none. See <see cref="MslFragmentIon.NeutralLoss"/>.
	/// </summary>
	public double NeutralLoss { get; }

	/// <summary>
	/// True when the fragment should be excluded from quantification.
	/// See <see cref="MslFragmentIon.ExcludeFromQuant"/>.
	/// </summary>
	public bool ExcludeFromQuant { get; }

	/// <summary>True for an internal ion (<see cref="SecondaryProductType"/> is not null).</summary>
	public bool IsInternalFragment => SecondaryProductType != null;

	/// <summary>True for a diagnostic ion (<see cref="ProductType"/> is <see cref="ProductType.D"/>).</summary>
	public bool IsDiagnosticIon => ProductType == ProductType.D;
}
