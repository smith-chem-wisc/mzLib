using UsefulProteomicsDatabases;

namespace Readers
{
    /// <summary>What kind of file a MassIVE run is.</summary>
    public enum MassiveRunKind
    {
        /// <summary>A vendor acquisition: a .raw, .wiff, .wiff2 or .d (file, folder or zipped).</summary>
        Vendor,

        /// <summary>A .zip in the raw collection, when the deposit has no vendor file (MSV000098930's 508).</summary>
        Archive,

        /// <summary>The submitter's own peak list (.mzML, .mzXML, .mgf, or one of them gzipped), when there is nothing else.</summary>
        PeakList,
    }

    /// <summary>
    /// One run of a MassIVE deposit: the file an SDRF row names, with the files that belong to the same acquisition.
    /// </summary>
    public sealed class MassiveRun
    {
        internal MassiveRun(string path, MassiveRunKind kind, bool isFolder, IReadOnlyList<MassiveFile> files)
        {
            Path = path;
            Kind = kind;
            IsFolder = isFolder;
            Files = files;
        }

        /// <summary>
        /// The run's path below the dataset: its primary file ("raw/Data/x.wiff"), or, for a folder run, the folder
        /// ("raw/x.d"), which a Bruker or Waters acquisition is.
        /// </summary>
        public string Path { get; }

        /// <summary>The last segment of <see cref="Path"/>: the name an SDRF's <c>comment[data file]</c> gives the run.</summary>
        public string FileName => Path[(Path.LastIndexOf('/') + 1)..];

        /// <summary>What kind of file the run is.</summary>
        public MassiveRunKind Kind { get; }

        /// <summary>Whether the run is a folder (a .d or a Waters .raw), listed file by file.</summary>
        public bool IsFolder { get; }

        /// <summary>The primary file, or null for a folder run.</summary>
        public MassiveFile? Primary => IsFolder ? null : Files[0];

        /// <summary>
        /// Every file of the run: the primary first, then its members (sidecars such as .wiff.scan and
        /// .timeseries.data, and other formats of the same acquisition), in listing order. For a folder run, every
        /// file below the folder.
        /// </summary>
        public IReadOnlyList<MassiveFile> Files { get; }

        /// <summary>The files beside the primary; for a folder run, all of <see cref="Files"/>.</summary>
        public IReadOnlyList<MassiveFile> Members => IsFolder ? Files : Files.Skip(1).ToList();

        /// <summary>The MassIVE collection the run is in, as <see cref="MassiveFile.Collection"/>.</summary>
        public string? Collection => Files[0].Collection;

        /// <summary>The update the run came in, or null, as <see cref="MassiveFile.UpdateName"/>.</summary>
        public string? UpdateName => Files[0].UpdateName;

        /// <summary>The listed size of all of <see cref="Files"/>.</summary>
        public long SizeBytes => Files.Sum(f => f.SizeBytes ?? 0);

        /// <inheritdoc />
        public override string ToString() => Path;
    }

    /// <summary>
    /// Which files of a MassIVE listing are the deposit's runs: the rule sdrf's drafter (MV2) and the aging pipeline
    /// read a MassIVE deposit by (<c>design/M1.md</c> section 7, agreed in sdrf 033, <c>ccms_peak</c> dropped 2026-11-10).
    /// </summary>
    /// <remarks>
    /// <list type="number">
    /// <item><b>Folder runs.</b> A folder whose name ends ".d" (Bruker, Agilent) or ".raw" (Waters) is ONE run, with
    /// every file below it.</item>
    /// <item><b>Stems.</b> Other files group into one acquisition by folder and stem. The stem is
    /// <see cref="SdrfFileNamePattern.Stem"/>'s, so a run is named exactly as the SDRF drafter names its files, after
    /// ".timeseries.data" (a SCIEX OS sidecar the drafter's list does not hold) is stripped.</item>
    /// <item><b>Vendor acquisitions.</b> A group holding a .raw, .wiff, .wiff2, .d, .raw.zip or .d.zip is a run, its
    /// primary the .wiff over the .wiff2 over the .raw (the curated SDRF of PXD038339 names the .wiff on all 26 rows).
    /// A .wiff.scan or .timeseries.data is a sidecar, never a run.</item>
    /// <item><b>Which are the runs.</b> Any vendor acquisition: those, and nothing else. Else any other .zip in the
    /// raw collection: each zip. Else the submitter's peak lists. Else none, which is a real answer (a search-only
    /// reanalysis), not a failure.</item>
    /// </list>
    /// <c>ccms_peak</c>, MassIVE's own converted copies, is never a run: MassIVE's FTP server will not open it, so it
    /// is only a link in the listing. A deposit whose only spectra are there (MSV000080763) has no runs.
    /// </remarks>
    public static class MassiveRuns
    {
        private const string TimeseriesSidecar = ".timeseries.data";

        // Vendor files, in the order the primary is chosen in when one acquisition has several.
        private static readonly string[] VendorExtensions = { ".wiff", ".wiff2", ".raw", ".d", ".raw.zip", ".d.zip" };

        // Peak lists, in the order the primary is chosen in.
        private static readonly string[] PeakListExtensions =
            { ".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz" };

        private static readonly string[] FolderRunExtensions = { ".d", ".raw" };

        /// <summary>The runs of a deposit, in listing order. Links are ignored.</summary>
        /// <param name="listing">The deposit's listing, as <see cref="MassiveClient.GetDatasetFilesAsync"/> gives it.</param>
        /// <returns>The runs; empty when the deposit has none. Never null.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="listing"/> or an entry in it is null.</exception>
        public static IReadOnlyList<MassiveRun> Of(IEnumerable<MassiveFile> listing)
        {
            ArgumentNullException.ThrowIfNull(listing);
            var files = new List<MassiveFile>();
            foreach (MassiveFile? file in listing)
            {
                if (file == null)
                    throw new ArgumentNullException(nameof(listing), "The listing holds a null entry.");
                if (!file.IsLink)
                    files.Add(file);
            }

            // Listing positions, for listing order (a deposit can list 30,000 files: MSV000094486).
            var position = new Dictionary<MassiveFile, int>(ReferenceEqualityComparer.Instance);
            for (int i = 0; i < files.Count; i++)
                position[files[i]] = i;

            // 1. Folder runs, keyed by the folder's path.
            var folderRuns = new Dictionary<string, List<MassiveFile>>(StringComparer.OrdinalIgnoreCase);
            var folderOrder = new List<string>();
            // 2. Everything else, by folder and stem, in first-seen order.
            var groups = new Dictionary<string, List<MassiveFile>>(StringComparer.OrdinalIgnoreCase);
            var groupOrder = new List<string>();

            foreach (MassiveFile file in files)
            {
                string? folder = RunFolderOf(file.Path);
                if (folder != null)
                {
                    if (!folderRuns.TryGetValue(folder, out List<MassiveFile>? inFolder))
                    {
                        folderRuns[folder] = inFolder = new List<MassiveFile>();
                        folderOrder.Add(folder);
                    }
                    inFolder.Add(file);
                    continue;
                }

                string key = ParentOf(file.Path) + "\u0000" + StemOf(file.FileName);
                if (!groups.TryGetValue(key, out List<MassiveFile>? group))
                {
                    groups[key] = group = new List<MassiveFile>();
                    groupOrder.Add(key);
                }
                group.Add(file);
            }

            // 3. Vendor acquisitions: the folder runs, and each group with a vendor file.
            var vendor = folderOrder
                .Select(folder => (Position: position[folderRuns[folder][0]],
                                   Run: new MassiveRun(folder, MassiveRunKind.Vendor, isFolder: true, folderRuns[folder])))
                .ToList();
            foreach (string key in groupOrder)
            {
                MassiveFile? primary = PrimaryOf(groups[key], VendorExtensions);
                if (primary != null)
                    vendor.Add((position[primary], RunOf(primary, groups[key], MassiveRunKind.Vendor)));
            }
            if (vendor.Count > 0)
                return vendor.OrderBy(v => v.Position).Select(v => v.Run).ToList();

            // 4a. Else each .zip in the raw collection.
            List<MassiveRun> archives = files
                .Where(f => string.Equals(f.Collection, "raw", StringComparison.OrdinalIgnoreCase)
                            && f.FileName.EndsWith(".zip", StringComparison.OrdinalIgnoreCase))
                .Select(f => new MassiveRun(f.Path, MassiveRunKind.Archive, isFolder: false, new[] { f }))
                .ToList();
            if (archives.Count > 0)
                return archives;

            // 4b. Else the submitter's peak lists, one run per acquisition.
            var peakLists = new List<(int Position, MassiveRun Run)>();
            foreach (string key in groupOrder)
            {
                MassiveFile? primary = PrimaryOf(groups[key], PeakListExtensions);
                if (primary != null)
                    peakLists.Add((position[primary], RunOf(primary, groups[key], MassiveRunKind.PeakList)));
            }
            return peakLists.OrderBy(p => p.Position).Select(p => p.Run).ToList();
        }

        /// <summary>
        /// The stem a run's files share: <see cref="SdrfFileNamePattern.Stem"/>'s, after ".timeseries.data" is stripped,
        /// so "X.wiff", "X.wiff2", "X.wiff.scan" and "X.timeseries.data" are all "X".
        /// </summary>
        internal static string StemOf(string fileName)
        {
            if (fileName.Length > TimeseriesSidecar.Length && fileName.EndsWith(TimeseriesSidecar, StringComparison.OrdinalIgnoreCase))
                fileName = fileName[..^TimeseriesSidecar.Length];
            return SdrfFileNamePattern.Stem(fileName);
        }

        /// <summary>The outermost folder of <paramref name="path"/> that is a run (".d", ".raw"), or null.</summary>
        private static string? RunFolderOf(string path)
        {
            string[] segments = path.Split('/');
            for (int i = 0; i < segments.Length - 1; i++)
            {
                string segment = segments[i];
                if (FolderRunExtensions.Any(e => segment.Length > e.Length && segment.EndsWith(e, StringComparison.OrdinalIgnoreCase)))
                    return string.Join('/', segments.Take(i + 1));
            }
            return null;
        }

        private static string ParentOf(string path)
        {
            int slash = path.LastIndexOf('/');
            return slash < 0 ? string.Empty : path[..slash];
        }

        /// <summary>
        /// The file of <paramref name="group"/> with the first of <paramref name="extensions"/> it has, or null. A sidecar
        /// (".wiff.scan", ".timeseries.data") is never chosen: no vendor or peak-list extension is a suffix of either.
        /// </summary>
        private static MassiveFile? PrimaryOf(List<MassiveFile> group, string[] extensions)
        {
            foreach (string extension in extensions)
            {
                MassiveFile? match = group.FirstOrDefault(f =>
                    f.FileName.Length > extension.Length && f.FileName.EndsWith(extension, StringComparison.OrdinalIgnoreCase));
                if (match != null)
                    return match;
            }
            return null;
        }

        private static MassiveRun RunOf(MassiveFile primary, List<MassiveFile> group, MassiveRunKind kind) =>
            new(primary.Path, kind, isFolder: false, new[] { primary }.Concat(group.Where(f => !ReferenceEquals(f, primary))).ToList());
    }
}
