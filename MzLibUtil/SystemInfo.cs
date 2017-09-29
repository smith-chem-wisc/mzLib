using Microsoft.Win32;
using System;
using System.Diagnostics;
using System.Management;
using System.Text;

namespace MzLibUtil
{
    public static class SystemInfo
    {
        #region Public Methods

        public static string CompleteSystemInfo()
        {
            StringBuilder fullSystemString = new StringBuilder();

            fullSystemString.Append(SystemProse() + "\n");

            fullSystemString.Append(DotNet());
            fullSystemString.Append(MsFileReader_FileIo());
            fullSystemString.Append(MsFileReader_Fregistry());
            fullSystemString.Append(MsFileReader_XRawfile2());

            return fullSystemString.ToString();
        }

        public static string SystemProse()
        {
            StringBuilder fullSystemProse = new StringBuilder();

            fullSystemProse.Append("Data files were processed on a " + GetManufacturer() + " computer ");
            fullSystemProse.Append("using " + WindowsOperatingSystemVersion());
            fullSystemProse.Append(" with a " + GetCpuRegister());
            fullSystemProse.Append(" and " + ProcessorCount() + " cores ");
            fullSystemProse.Append("operating at " + GetMaxClockSpeed() + "GHz ");
            fullSystemProse.Append("and " + InstalledRam() + "GB installed RAM.");

            return fullSystemProse.ToString();
        }

        #endregion Public Methods

        #region Private Methods

        private static string GetManufacturer()
        {
            string computerModel = "UNDETERMINED";
            try
            {
                System.Management.SelectQuery query = new System.Management.SelectQuery(@"Select * from Win32_ComputerSystem");

                using (System.Management.ManagementObjectSearcher searcher = new System.Management.ManagementObjectSearcher(query))
                {
                    foreach (System.Management.ManagementObject process in searcher.Get())
                    {
                        process.Get();
                        computerModel = process["Manufacturer"].ToString() + " " + process["Model"].ToString();
                    }
                }
                return computerModel;
            }
            catch
            {
                return computerModel;
            }
        }

        private static string GetWindowsOs()
        {
            try
            {
                var reg = Registry.LocalMachine.OpenSubKey(@"SOFTWARE\Microsoft\Windows NT\CurrentVersion");

                return ((string)reg.GetValue("ProductName"));
            }
            catch
            {
                return "UNDETERMINED OPERATING SYSTEM";
            }
        }

        private static string GetCpuRegister()
        {
            try
            {
                if (Environment.Is64BitOperatingSystem)
                    return "64-Bit processor";
                else
                    return "32-Bit processor";
            }
            catch
            {
                return "UNDETERMINED-BIT PROCESSOR";
            }
        }

        private static string GetMaxClockSpeed()
        {
            try
            {
                RegistryKey registrykeyHKLM = Registry.LocalMachine;
                string keyPath = @"HARDWARE\DESCRIPTION\System\CentralProcessor\0";
                RegistryKey registrykeyCPU = registrykeyHKLM.OpenSubKey(keyPath, false);
                string MHz = registrykeyCPU.GetValue("~MHz").ToString();
                double numericalMHz = Convert.ToDouble(MHz) / 1000d;
                registrykeyCPU.Close();
                return numericalMHz.ToString();
            }
            catch
            {
                return "UNDETERMINED";
            }
        }

        private static string WindowsOperatingSystemVersion()
        {
            try
            {
                return ("Windows version " + Environment.OSVersion.Version.ToString());
            }
            catch
            {
                return "UNDETERMINED OPERATING SYSTEM";
            }
        }

        private static string DotNet()
        {
            try
            {
                const string subkey = @"SOFTWARE\Microsoft\NET Framework Setup\NDP\v4\Full\";

                using (RegistryKey ndpKey = RegistryKey.OpenBaseKey(RegistryHive.LocalMachine, RegistryView.Registry32).OpenSubKey(subkey))
                {
                    if (ndpKey != null && ndpKey.GetValue("Release") != null)
                    {
                        return ".NET Framework Version: " + CheckFor45PlusVersion((int)ndpKey.GetValue("Release"));
                    }
                    else
                    {
                        return ".NET Framework Version 4.5 or later is not detected.";
                    }
                }
            }
            catch
            {
                return "Windows .Net Version could not be determined.\n";
            }
        }

        private static string InstalledRam()
        {
            try
            {
                string Query = "SELECT Capacity FROM Win32_PhysicalMemory";
                ManagementObjectSearcher searcher = new ManagementObjectSearcher(Query);

                UInt64 Capacity = 0;
                foreach (ManagementObject WniPART in searcher.Get())
                {
                    Capacity += Convert.ToUInt64(WniPART.Properties["Capacity"].Value);
                }

                return ((Capacity / 1073741824).ToString());
            }
            catch
            {
                return "UNKNOWN ";
            }
        }

        private static string ProcessorCount()
        {
            try
            {
                return (Environment.ProcessorCount.ToString());
            }
            catch
            {
                return "AN UNKNOWN NUMBER OF ";
            }
        }

        private static string MsFileReader_FileIo()
        {
            try
            {
                FileVersionInfo myFileVersionInfo = FileVersionInfo.GetVersionInfo(@"C:\Program Files\Thermo\MSFileReader\Fileio_x64.dll");
                return ("Thermo MSFileReader " + myFileVersionInfo.FileDescription + "\t\t" +
                              "Version: " + myFileVersionInfo.FileVersion + "\n");
            }
            catch
            {
                return @"C:\Program Files\Thermo\MSFileReader\Fileio_x64.dll MISSING\n";
            }
        }

        private static string MsFileReader_Fregistry()
        {
            try
            {
                FileVersionInfo myFileVersionInfo = FileVersionInfo.GetVersionInfo(@"C:\Program Files\Thermo\MSFileReader\fregistry_x64.dll");
                return ("Thermo MSFileReader " + myFileVersionInfo.FileDescription + "\t\t" +
                              "Version: " + myFileVersionInfo.FileVersion + "\n");
            }
            catch
            {
                return @"C:\Program Files\Thermo\MSFileReader\fregistry_x64.dll MISSING\n";
            }
        }

        private static string MsFileReader_XRawfile2()
        {
            try
            {
                FileVersionInfo myFileVersionInfo = FileVersionInfo.GetVersionInfo(@"C:\Program Files\Thermo\MSFileReader\XRawfile2_x64.dll");
                return ("Thermo MSFileReader " + myFileVersionInfo.FileDescription + "\t\t" +
                              "Version: " + myFileVersionInfo.FileVersion + "\n");
            }
            catch
            {
                return @"C:\Program Files\Thermo\MSFileReader\XRawfile2_x64.dll MISSING\n";
            }
        }

        private static string CheckFor45PlusVersion(int releaseKey)
        {
            if (releaseKey >= 460798)
                return "4.7 or later";
            if (releaseKey >= 394802)
                return "4.6.2";
            if (releaseKey >= 394254)
            {
                return "4.6.1";
            }
            if (releaseKey >= 393295)
            {
                return "4.6";
            }
            if ((releaseKey >= 379893))
            {
                return "4.5.2";
            }
            if ((releaseKey >= 378675))
            {
                return "4.5.1";
            }
            if ((releaseKey >= 378389))
            {
                return "4.5";
            }
            // This code should never execute. A non-null release key should mean
            // that 4.5 or later is installed.
            return "No 4.5 or later version detected";
        }

        #endregion Private Methods
    }
}