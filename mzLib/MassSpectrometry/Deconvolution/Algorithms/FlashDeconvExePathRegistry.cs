using System;
using System.Collections.Concurrent;
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// Static registry of validated FLASHDeconv executable paths.
    ///
    /// Each <see cref="Deconvoluter.Deconvolute(MzSpectrum, DeconvolutionParameters, MzLibUtil.MzRange)"/>
    /// call constructs a fresh <see cref="RealFLASHDeconvolutionAlgorithm"/>, and
    /// resolving the FLASHDeconv exe path costs File.Exists syscalls (the explicit
    /// path, plus a walk of well-known paths and the PATH env var on misses).
    /// On a per-scan hot path that's a few-hundred-ms stack of redundant syscalls
    /// per minute. The exe doesn't move between calls, so first-call validation
    /// is enough -- this registry caches the resolved path so subsequent calls
    /// skip the filesystem.
    ///
    /// Production callers that know the path up front (typically MetaMorpheus via
    /// <c>GlobalSettings.FLASHDeconvExecutablePath</c>) should call
    /// <see cref="Register"/> once at startup. Ad-hoc callers don't have to do
    /// anything: the algorithm caches lazily on first resolution either way.
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
        /// Validate the given exe path once and cache it so subsequent
        /// deconvolution calls skip the filesystem check.
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
        /// Cache a resolution that has already been validated (typically by the
        /// algorithm's own resolve pass). Keyed identically to
        /// <see cref="TryGet"/>: explicit path or default-search sentinel.
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
