using System;
using System.Collections.Generic;

namespace Readers.SpectralLibrary;

/// <summary>
/// One confident identification from one raw file, as <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>
/// reads it.
/// </summary>
/// <param name="FullSequence">
///   The modified sequence, in the notation the library's entries use (it is matched to an entry
///   by exact <see cref="Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry.FullSequence"/>
///   and charge).
/// </param>
/// <param name="BaseSequence">The unmodified sequence.</param>
/// <param name="Charge">The precursor charge.</param>
/// <param name="RetentionTimeMinutes">The scan's retention time in the run, in minutes.</param>
/// <remarks>
/// The caller passes only confident target identifications: no decoys or contaminants, a q-value of
/// at most 1%, one unambiguous sequence, and no match-between-runs detections (their time was placed
/// by alignment, not measured). A time that is not finite, or is negative, is ignored.
/// </remarks>
public sealed record MslRunPsm(string FullSequence, string BaseSequence, int Charge, double RetentionTimeMinutes);

/// <summary>
/// The confident identifications of one raw file (one LC-MS run), for
/// <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>.
/// </summary>
/// <param name="RunName">The raw file's name, used in the report.</param>
/// <param name="Psms">The run's confident identifications.</param>
/// <param name="RunLengthMinutes">
///   The run's length in minutes, which sets the void window at its start. When unknown (NaN), the
///   latest identification time stands in for it.
/// </param>
public sealed record MslRunObservations(string RunName, IReadOnlyList<MslRunPsm> Psms, double RunLengthMinutes = double.NaN);

/// <summary>
/// Thresholds for <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>. The defaults are the values
/// measured and tested on PXReprise's rat, mouse and human library chains (2026-11).
/// </summary>
/// <param name="AnchorsPerRun">Calibration anchors chosen per run (fewer when the run has fewer candidates).</param>
/// <param name="MinimumAnchors">A run calibrated on fewer anchors than this contributes no times.</param>
/// <param name="MaxResidualSd">A run whose calibration residual SD (iRT) is larger contributes no times.</param>
/// <param name="AgreementTolerance">Two times within this many iRT agree.</param>
/// <param name="StrongDisagreement">Two times more than this many iRT apart strongly disagree. Between the two tolerances the result is inconclusive.</param>
/// <param name="PredictionTolerance">How far (iRT) a time seen in a single run may be from the length-corrected prediction.</param>
/// <param name="LongPeptidePredictionTolerance">The same, for unmodified peptides over 30 residues checked against the long-peptide predictor.</param>
/// <param name="DatasetDisagreementLimit">
///   When more than this share of a dataset's times that overlap measured library times disagree with
///   them, the dataset only confirms: it fills no time and challenges none.
/// </param>
/// <param name="VoidFraction">Times in this first fraction of the run (the void) are refused.</param>
/// <param name="DoublePeakGap">A precursor whose identifications in a run have a gap larger than this (iRT) is a double peak, and gives no time from that run.</param>
/// <param name="StreakSpan">A precursor whose identifications in a run span more than this (iRT) is a streak, and gives no time from that run.</param>
public sealed record MslRetentionTimeUpdateOptions(
	int AnchorsPerRun = 1000,
	int MinimumAnchors = 100,
	double MaxResidualSd = 10.0,
	double AgreementTolerance = 7.0,
	double StrongDisagreement = 12.0,
	double PredictionTolerance = 15.0,
	double LongPeptidePredictionTolerance = 22.0,
	double DatasetDisagreementLimit = 0.40,
	double VoidFraction = 0.04,
	double DoublePeakGap = 3.0,
	double StreakSpan = 5.0);

/// <summary>What one run contributed to a retention-time update.</summary>
/// <param name="RunName">The raw file's name.</param>
/// <param name="AnchorsSelected">Anchors chosen for the run.</param>
/// <param name="AnchorsPredicted">Chosen anchors that had a prediction (the ones fitted).</param>
/// <param name="Gate">Why the run contributed no times ("too few anchors", "fit failed", "residual SD"), or null.</param>
/// <param name="ResidualSd">The calibration's residual SD in iRT; NaN when not fitted.</param>
/// <param name="TimesAccepted">Precursors whose time from this run was accepted.</param>
/// <param name="RefusedDoublePeak">Refused: identifications in two separated places.</param>
/// <param name="RefusedStreak">Refused: identifications spread over a long stretch.</param>
/// <param name="RefusedVoid">Refused: in the void at the start of the run.</param>
/// <param name="RefusedExtrapolated">Refused: outside the calibrated range.</param>
/// <param name="RefusedUnreliableStretch">Refused: in a third of the run the calibration does not fit well.</param>
public sealed record MslRunCalibrationReport(
	string RunName,
	int AnchorsSelected,
	int AnchorsPredicted,
	string? Gate,
	double ResidualSd,
	int TimesAccepted,
	int RefusedDoublePeak,
	int RefusedStreak,
	int RefusedVoid,
	int RefusedExtrapolated,
	int RefusedUnreliableStretch)
{
	/// <summary>True when the run was calibrated and could contribute times.</summary>
	public bool Calibrated => Gate is null;
}

/// <summary>The outcome of <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>.</summary>
public sealed record MslRetentionTimeUpdateResult
{
	/// <summary>The written library.</summary>
	public required string OutputPath { get; init; }

	/// <summary>One report per run, in input order.</summary>
	public required IReadOnlyList<MslRunCalibrationReport> Runs { get; init; }

	/// <summary>Missing times filled with the dataset's time (status observed once).</summary>
	public int Filled { get; init; }

	/// <summary>
	/// Missing times filled from a single run whose time the prediction disagreed with: stored as
	/// measured, but challenged and struck by the prediction.
	/// </summary>
	public int FilledChallenged { get; init; }

	/// <summary>Times newly confirmed (the time becomes the mean of the two).</summary>
	public int Confirmed { get; init; }

	/// <summary>Times already confirmed that agreed again (unchanged).</summary>
	public int AlreadyConfirmed { get; init; }

	/// <summary>Times challenged (kept, but marked) by a strong disagreement.</summary>
	public int Challenged { get; init; }

	/// <summary>Challenged times replaced by a second strong disagreement (status observed once).</summary>
	public int ReplacedSecondStrike { get; init; }

	/// <summary>Of those, replacements of a prediction-struck time that the prediction challenged again.</summary>
	public int ReplacedAndStruckAgain { get; init; }

	/// <summary>
	/// Times replaced by the prediction tiebreak: the stored time was more than 15 iRT from the
	/// prediction and the new one within 5.
	/// </summary>
	public int ReplacedTiebreak { get; init; }

	/// <summary>Times between the two tolerances: no change.</summary>
	public int Inconclusive { get; init; }

	/// <summary>Single-run times for an entry with a stored time that the prediction refused.</summary>
	public int RefusedSingleRunTimes { get; init; }

	/// <summary>Dataset times held back because the dataset only confirms (see <see cref="OnlyConfirmations"/>).</summary>
	public int HeldBack { get; init; }

	/// <summary>Dataset times for precursors (full sequence and charge) not in the library.</summary>
	public int NotInLibrary { get; init; }

	/// <summary>Dataset times that overlap a measured library time.</summary>
	public int Overlaps { get; init; }

	/// <summary>Of those, how many disagree with the stored time.</summary>
	public int OverlapsDisagreeing { get; init; }

	/// <summary>
	/// True when too many overlaps disagreed (<see cref="MslRetentionTimeUpdateOptions.DatasetDisagreementLimit"/>),
	/// so the dataset only confirmed.
	/// </summary>
	public bool OnlyConfirmations { get; init; }

	/// <summary>Distinct sequences sent to the iRT predictor.</summary>
	public int PredictionsRequested { get; init; }

	/// <summary>Of those, how many came back without a value.</summary>
	public int PredictionsMissing { get; init; }

	/// <summary>Library entries whose time or status changed.</summary>
	public int EntriesChanged { get; init; }
}
