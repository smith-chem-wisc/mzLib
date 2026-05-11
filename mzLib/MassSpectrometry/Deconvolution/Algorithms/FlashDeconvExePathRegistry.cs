#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// Caller-side helper for locating and caching the FLASHDeconv executable.
    ///
    /// <see cref="RealFLASHDeconvolutionAlgorithm"/> deliberately does NOT search
    /// for the executable; it uses whatever path is on its parameters object and
    /// throws if that path is missing or invalid. This class is where the
    /// "go find FLASHDeconv on this machine" logic lives, so that exe discovery
    /// is an explicit step the caller performs once -- not a side effect of
    /// running deconvolution.
    ///
    /// Typical usage (e.g. MetaMorpheus at startup):
    /// <code>
    /// string exe = FlashDeconvExePathRegistry.Resolve(GlobalSettings.FLASHDeconvExecutablePath);
    /// // ...later, when configuring deconvolution params...
    /// var p = new RealFLASHDeconvolutionParameters(flashDeconvExePath: exe);
    /// </code>
    ///
    /// The registry caches resolved paths so repeated <see cref="Resolve"/>
    /// calls skip the filesystem walk.
    /// </summary>
    public static class FlashDeconvExePathRegistry
    {
        // Sentinel key for "no explicit path was supplied; this is the result of
        // walking the well-known list + PATH". A real path string can never collide
        // with this because the inner '<' / '>' aren't legal in filesystem paths.
        internal const string DefaultSearchSentinel = "<default-search>";

        private static readonly ConcurrentDictionary<string, string> _validated
            = new ConcurrentDictionary<string, string>();

        /// <summary>
        /// Hardcoded install paths probed by <see cref="Resolve(string?)"/> when no
        /// explicit FLASHDeconv path is supplied and PATH search misses.
        /// </summary>
        public static readonly IReadOnlyList<string> DefaultWellKnownPaths = new[]
        {
            @"C:\Program Files\OpenMS-3.0.0-pre-HEAD-2023-06-17\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.5.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.4.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.3.0\bin\FLASHDeconv.exe",
            @"C:\Program Files\OpenMS-3.0.0\bin\FLASHDeconv.exe",
            "/usr/bin/FLASHDeconv",
            "/usr/local/bin/FLASHDeconv",
            "/opt/openms/bin/FLASHDeconv",
        };

        /// <summary>
        /// Register an already-known FLASHDeconv path. Validates that the file
        /// exists, then caches the result so subsequent <see cref="Resolve(string?)"/>
        /// calls with the same path return immediately.
        /// </summary>
        /// <exception cref="ArgumentException">path is null/whitespace.</exception>
        /// <exception cref="FileNotFoundException">path does not exist on disk.</exception>
        public static void Register(string path)
        {
            if (string.IsNullOrWhiteSpace(path))
                throw new ArgumentException("Path must be a non-empty string.", nameof(path));
            if (!File.Exists(path))
                throw new FileNotFoundException(
                    $"FLASHDeconv not found at: {path}", path);
            _validated[path] = path;
        }

        /// <summary>
        /// Locate the FLASHDeconv executable. Search order:
        ///   1. <paramref name="explicitPath"/> if supplied
        ///   2. <see cref="DefaultWellKnownPaths"/>
        ///   3. directories on the PATH environment variable
        /// Caches the result so repeated calls skip the filesystem walk.
        /// </summary>
        /// <param name="explicitPath">Optional explicit path; pass null/empty to
        /// trigger the default search.</param>
        /// <returns>The validated absolute path of FLASHDeconv.</returns>
        /// <exception cref="FileNotFoundException">
        /// <paramref name="explicitPath"/> was supplied but does not exist, or no
        /// explicit path was supplied and the default search found nothing.
        /// </exception>
        public static string Resolve(string? explicitPath = null)
        {
            if (TryGet(explicitPath, out string cached))
                return cached;

            string resolved = Resolve(
                explicitPath,
                DefaultWellKnownPaths,
                Environment.GetEnvironmentVariable("PATH"));

            CacheValidated(explicitPath, resolved);
            return resolved;
        }

        /// <summary>
        /// Test-friendly overload that takes the well-known-paths list and PATH-env
        /// value as parameters instead of reading them from compiled-in defaults +
        /// the live process environment. Lets tests deterministically exercise the
        /// well-known-search, PATH-search, and not-found-anywhere branches without
        /// depending on what's installed on the runner.
        /// Does NOT consult or update the cache; pure resolution.
        /// </summary>
        internal static string Resolve(
            string? explicitPath,
            IEnumerable<string> wellKnownPaths,
            string? pathEnv)
        {
            if (!string.IsNullOrWhiteSpace(explicitPath))
            {
                if (File.Exists(explicitPath)) return explicitPath!;
                throw new FileNotFoundException(
                    $"FLASHDeconv not found at: {explicitPath}", explicitPath);
            }

            foreach (string c in wellKnownPaths)
                if (File.Exists(c)) return c;

            if (pathEnv != null)
            {
                string[] names = OperatingSystem.IsWindows()
                    ? new[] { "FLASHDeconv.exe" }
                    : new[] { "FLASHDeconv", "flashdeconv" };
                foreach (string dir in pathEnv.Split(Path.PathSeparator))
                    foreach (string name in names)
                    {
                        string full = Path.Combine(dir, name);
                        if (File.Exists(full)) return full;
                    }
            }

            throw new FileNotFoundException(
                "FLASHDeconv not found. Set RealFLASHDeconvolutionParameters.FLASHDeconvExePath " +
                "or install OpenMS to a well-known location.");
        }

        /// <summary>
        /// Number of entries currently cached. Useful for tests that want to
        /// assert registration happened without depending on internal state.
        /// </summary>
        public static int Count => _validated.Count;

        /// <summary>
        /// Forget every cached path. Test-only -- production code has no reason
        /// to clear because validation is permissive (a stale cached entry just
        /// pushes the actual existence check to Process.Start, which surfaces
        /// the same failure with a clear error).
        /// </summary>
        internal static void Clear() => _validated.Clear();

        /// <summary>
        /// Look up a previously-validated resolution. Key is the explicit path
        /// string, or the default-search sentinel when <paramref name="explicitPath"/>
        /// is null/whitespace (i.e. the caller wants the algorithm to search
        /// well-known paths + PATH itself).
        /// </summary>
        internal static bool TryGet(string? explicitPath, out string resolved)
        {
            string key = string.IsNullOrWhiteSpace(explicitPath)
                ? DefaultSearchSentinel
                : explicitPath!;
            return _validated.TryGetValue(key, out resolved!);
        }

        /// <summary>
        /// Cache a resolution that has already been validated. Keyed identically
        /// to <see cref="TryGet"/>: explicit path or default-search sentinel.
        /// </summary>
        internal static void CacheValidated(string? explicitPath, string resolved)
        {
            string key = string.IsNullOrWhiteSpace(explicitPath)
                ? DefaultSearchSentinel
                : explicitPath!;
            _validated[key] = resolved;
        }
    }
}
