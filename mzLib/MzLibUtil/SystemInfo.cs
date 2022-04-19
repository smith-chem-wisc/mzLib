using Easy.Common;
using Microsoft.Win32;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.Versioning;
using System.Text;

namespace MzLibUtil
{
    public static class SystemInfo
    {
        public static string CompleteSystemInfo()
        {
            StringBuilder fullSystemString = new StringBuilder();

            fullSystemString.Append(SystemProse() + "\n");

            fullSystemString.Append(".NET version: " + GetDotNetVersion());

            return fullSystemString.ToString();
        }

        public static string SystemProse()
        {
            StringBuilder fullSystemProse = new StringBuilder();

            fullSystemProse.Append("Data files were processed on a computer ");
            fullSystemProse.Append("running " + GetOperatingSystem());
            fullSystemProse.Append(" with a " + GetCpuRegister());
            fullSystemProse.Append(" " + GetProcessorName() + " processor");
            fullSystemProse.Append(" with " + ProcessorCount() + " threads ");
            
            fullSystemProse.Append("and " + InstalledRam() + "GB installed RAM.");

            return fullSystemProse.ToString();
        }

        private static string GetOperatingSystem()
        {
            try
            {
                return RuntimeInformation.OSDescription;
            }
            catch
            {
                return "UNKNOWN OPERATING SYSTEM";
            }
        }

        private static string GetCpuRegister()
        {
            try
            {
                if (Environment.Is64BitOperatingSystem)
                    return "64-bit";
                else
                    return "32-bit";
            }
            catch
            {
                return "UNDETERMINED-BIT";
            }
        }

        /// <summary>
        /// Gets the name of the machine's processor
        /// // see https://github.com/NimaAra/Easy.Common/blob/master/Easy.Common/DiagnosticReport/DiagnosticReport.cs
        /// </summary>
        private static string GetProcessorName()
        {
            string processorName = "UNKNOWN";

            try
            {
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                {
                    var key = Registry.LocalMachine.OpenSubKey(@"HARDWARE\DESCRIPTION\System\CentralProcessor\0\");
                    processorName = key?.GetValue("ProcessorNameString").ToString() ?? "Not Found";
                }

                if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                {
                    const string CPUFile = "/proc/cpuinfo";
                    var cpuLine = File.ReadLines(CPUFile)
                        .FirstOrDefault(l => l.StartsWith("model name", StringComparison.InvariantCultureIgnoreCase));

                    if (cpuLine != null)
                    {
                        const string Separator = ": ";
                        var startIdx = cpuLine.IndexOf(Separator, StringComparison.Ordinal) + Separator.Length;
                        processorName = cpuLine.Substring(startIdx, cpuLine.Length - startIdx);
                    }
                }

                if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                {
                    processorName = AsBashCommand("sysctl -n machdep.cpu.brand_string").TrimEnd();
                }
            }
            catch
            {

            }

            return processorName;
        }

        /// <summary>
        /// see https://github.com/NimaAra/Easy.Common/blob/master/Easy.Common/DiagnosticReport/DiagnosticReport.cs
        /// </summary>
        private static string AsBashCommand(string command)
        {
            var escapedArgs = command.Replace("\"", "\\\"");
            var p = new Process
            {
                EnableRaisingEvents = true,
                StartInfo = new ProcessStartInfo
                {
                    FileName = "sh",
                    Arguments = $"-c \"{escapedArgs}\"",
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                    UseShellExecute = false,
                    CreateNoWindow = true
                }
            };

            p.Start();
            var result = p.StandardOutput.ReadToEnd();
            p.WaitForExit();
            return result;
        }

        private static string GetDotNetVersion()
        {
            try
            {
                return RuntimeInformation.FrameworkDescription;
            }
            catch
            {
                return "UNKNOWN .NET VERSION";
            }
        }

        /// <summary>
        /// see https://github.com/NimaAra/Easy.Common/blob/master/Easy.Common/DiagnosticReport/DiagnosticReport.cs
        /// </summary>
        private static string InstalledRam()
        {
            string amountOfRamInGb = "UNKNOWN ";

            try
            {
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                {
                    GetPhysicallyInstalledSystemMemory(out var installedMemoryKb);
                    amountOfRamInGb = ((long)UnitConverter.KiloBytesToMegaBytes(installedMemoryKb).MegaBytesToGigaBytes()).ToString("F0");
                }

                if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                {
                    const string MemFile = "/proc/meminfo";
                    var memLine = File.ReadLines(MemFile)
                        .FirstOrDefault(l => l.StartsWith("MemTotal:", StringComparison.InvariantCultureIgnoreCase));

                    if (memLine != null)
                    {
                        const string BeginSeparator = ":";
                        const string EndSeparator = "kB";
                        var startIdx = memLine.IndexOf(BeginSeparator, StringComparison.Ordinal) + BeginSeparator.Length;
                        var endIdx = memLine.IndexOf(EndSeparator, StringComparison.Ordinal);
                        var memStr = memLine.Substring(startIdx, endIdx - startIdx);
                        amountOfRamInGb = (long.Parse(memStr) / 1000_000).ToString();
                    }
                }

                if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                {
                    var memStr = AsBashCommand("sysctl -n hw.memsize");
                    amountOfRamInGb = (long.Parse(memStr) / 1000_000_000).ToString();
                }
            }
            catch
            {
                
            }

            return amountOfRamInGb;
        }

        /// <summary>
        /// see https://github.com/NimaAra/Easy.Common/blob/master/Easy.Common/DiagnosticReport/DiagnosticReport.cs
        /// <see href="https://msdn.microsoft.com/en-us/library/windows/desktop/cc300158(v=vs.85).aspx"/>
        /// </summary>
        [DllImport("kernel32.dll")]
        [return: MarshalAs(UnmanagedType.Bool)]
        private static extern bool GetPhysicallyInstalledSystemMemory(out long totalMemoryInKilobytes);

        private static string ProcessorCount()
        {
            try
            {
                return (Environment.ProcessorCount.ToString());
            }
            catch
            {
                return "AN UNKNOWN NUMBER OF";
            }
        }
    }
}