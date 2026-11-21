using System;
using System.IO;
using System.Reflection;

namespace UsefulProteomicsDatabases.Annotation
{
    /// <summary>The "#!key value" header lines the annotation tables (GO, member resolution, Reactome) share.</summary>
    internal static class TsvHeader
    {
        /// <summary>Writes "#!key value"; a null value omits the line.</summary>
        public static void Write(TextWriter output, string key, string value)
        {
            if (value == null)
            {
                return;
            }
            RejectSeparators(value, key);
            output.Write($"#!{key} {value}" + "\n");
        }

        /// <summary>
        /// Writes "#!mzlib_version" (this assembly's informational version, commit included) and
        /// "#!mzlib_release" (the released version, or "none" for a local build), so an ingester can refuse
        /// a file no release produced without parsing the version itself.
        /// </summary>
        public static void WriteProducer(TextWriter output)
        {
            string version = typeof(TsvHeader).Assembly
                .GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;
            Write(output, "mzlib_version", version);
            Write(output, "mzlib_release", ReleaseOf(version));
        }

        /// <summary>
        /// The release an informational version names, or "none". release.yml builds with
        /// /p:Version=(git tag); a build without it gets the SDK default 1.0.0, which no release carries.
        /// </summary>
        public static string ReleaseOf(string informationalVersion)
        {
            string version = informationalVersion?.Split('+')[0];
            return string.IsNullOrEmpty(version) || version == "1.0.0" ? "none" : version;
        }

        public static void RejectSeparators(string value, string field)
        {
            if (value != null && value.IndexOfAny(new[] { '\t', '\n', '\r' }) >= 0)
            {
                throw new ArgumentException($"The value of '{field}' contains a tab or line break, which would shift the row.");
            }
        }
    }
}
