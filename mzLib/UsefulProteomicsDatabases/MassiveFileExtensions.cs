using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Name rules over a MassIVE listing (<see cref="MassiveClient.GetDatasetFilesAsync"/>): which files are a deposited
    /// SDRF, which may be an experimental-design table, and which are Spectronaut output. Names only: nothing is
    /// downloaded or read.
    /// </summary>
    public static class MassiveFileExtensions
    {
        private static readonly string[] SdrfExtensions = { ".tsv", ".txt", ".xlsx" };

        /// <summary>
        /// Whether <paramref name="file"/> is an SDRF the submitter deposited: its name contains "sdrf" (any case) and
        /// ends .tsv, .txt or .xlsx, in any collection. A link is never one.
        /// </summary>
        /// <remarks>
        /// Checked on 2026-11-03 against sdrf's records (sdrf 026, 033): it finds all seven deposited SDRFs known, and
        /// matches none of the 49 design tables in sdrf's 121-deposit inventory, nor MSV000103504's
        /// <c>experiment_annotation.tsv</c> beside its <c>sdrf.tsv</c>. A ".csv" SDRF has not been seen, so it is not one.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="file"/> is null.</exception>
        public static bool IsDepositedSdrf(this MassiveFile file)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));
            return !file.IsLink
                && file.FileName.Contains("sdrf", StringComparison.OrdinalIgnoreCase)
                && SdrfExtensions.Any(e => file.FileName.EndsWith(e, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Whether <paramref name="file"/> may be an experimental-design table, by the names sdrf found them under
        /// (sdrf 028): <c>*ConditionSetup*.tsv</c> (Spectronaut), <c>*annotation*.csv|tsv|txt</c>,
        /// <c>experimentalDesign*.txt</c> and <c>experimental.design.txt</c>, any case, in any collection. A deposited
        /// SDRF is not one, and neither is a link.
        /// </summary>
        /// <remarks>
        /// Only a candidate: whether it IS a design table is decided by its header, by the reader that opens it. The
        /// tables are not only in <c>metadata/</c> (MSV000090777 keeps 27 copies under <c>quant/</c>,
        /// <c>quant_stats/</c> and <c>search/</c>), so every collection is looked at.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="file"/> is null.</exception>
        public static bool IsDesignTableCandidate(this MassiveFile file)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));
            if (file.IsLink || file.IsDepositedSdrf())
                return false;

            string name = file.FileName;
            bool Ends(string extension) => name.EndsWith(extension, StringComparison.OrdinalIgnoreCase);
            bool Has(string part) => name.Contains(part, StringComparison.OrdinalIgnoreCase);

            return (Has("ConditionSetup") && Ends(".tsv"))
                || (Has("annotation") && (Ends(".csv") || Ends(".tsv") || Ends(".txt")))
                || (name.StartsWith("experimentalDesign", StringComparison.OrdinalIgnoreCase) && Ends(".txt"))
                || string.Equals(name, "experimental.design.txt", StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>The deposited SDRFs in a listing, in listing order.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="files"/> is null.</exception>
        public static IEnumerable<MassiveFile> WhereDepositedSdrf(this IEnumerable<MassiveFile> files) =>
            (files ?? throw new ArgumentNullException(nameof(files))).Where(f => f.IsDepositedSdrf());

        /// <summary>The design-table candidates in a listing, in listing order.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="files"/> is null.</exception>
        public static IEnumerable<MassiveFile> WhereDesignTableCandidate(this IEnumerable<MassiveFile> files) =>
            (files ?? throw new ArgumentNullException(nameof(files))).Where(f => f.IsDesignTableCandidate());

        private static readonly string[] SpectronautExtensions = { ".sne", ".kit" };

        /// <summary>
        /// Whether <paramref name="file"/> is Spectronaut output: an analysis (<c>.sne</c>) or a spectral library
        /// (<c>.kit</c>), any case, in any collection. A link is never one.
        /// </summary>
        /// <remarks>
        /// A signal for telling DIA deposits from DDA ones by their files, because a record's text understates DIA
        /// (aging 037, 2026-11-10): MSV000087085, MSV000096075, MSV000098115 and MSV000094212 carry Spectronaut files while
        /// their records name no DIA. Spectronaut is DIA software, but a library is often built from DDA runs
        /// (MSV000090777's <c>.kit</c> is a "DDA_Zeno_120min_spectral_library"), so this says the deposit was analysed
        /// with Spectronaut, not that every run is DIA.
        /// <para>
        /// Any collection, not only <c>search/</c>, <c>quant/</c> and <c>library/</c>: MSV000096226 keeps one of its two
        /// <c>.sne</c> files under <c>raw/</c>.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="file"/> is null.</exception>
        public static bool IsSpectronautOutput(this MassiveFile file)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));
            return !file.IsLink
                && SpectronautExtensions.Any(e => file.FileName.Length > e.Length && file.FileName.EndsWith(e, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>The Spectronaut output in a listing, in listing order.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="files"/> is null.</exception>
        public static IEnumerable<MassiveFile> WhereSpectronautOutput(this IEnumerable<MassiveFile> files) =>
            (files ?? throw new ArgumentNullException(nameof(files))).Where(f => f.IsSpectronautOutput());

        /// <summary>
        /// Whether a deposit's listing holds any Spectronaut output: the deposit was analysed with Spectronaut, a DIA
        /// signal its record's text may not give. See <see cref="IsSpectronautOutput"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="files"/> is null.</exception>
        public static bool HasSpectronautOutput(this IEnumerable<MassiveFile> files) => files.WhereSpectronautOutput().Any();
    }
}
