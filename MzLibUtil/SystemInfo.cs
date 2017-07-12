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

            fullSystemString.Append(GetWindowsOs());
            fullSystemString.Append(GetCpuRegister());
            fullSystemString.Append(WindowsOperatingSystemVersion());
            fullSystemString.Append(DotNet());
            fullSystemString.Append(InstalledRam());
            fullSystemString.Append(ProcessorCount());
            fullSystemString.Append(MsFileReader_FileIo());
            fullSystemString.Append(MsFileReader_Fregistry());
            fullSystemString.Append(MsFileReader_XRawfile2());

            return fullSystemString.ToString();
        }

        #endregion Public Methods

        #region Private Methods
        
        private static string GetWindowsOs()
        {
            try
            {
                var reg = Registry.LocalMachine.OpenSubKey(@"SOFTWARE\Microsoft\Windows NT\CurrentVersion");

                return ((string)reg.GetValue("ProductName") +"\n");
            }
            catch
            {
                return "Windows operating system could not be determined.\n";
            }
        }

        private static string GetCpuRegister()
        {
            try
            {
                if (Environment.Is64BitOperatingSystem)
                    return "64-Bit\n";
                else
                    return "32-Bit\n";
            }
            catch
            {
                return "CPU register could not be determined.\n";
            }
        }

        private static string WindowsOperatingSystemVersion()
        {
            try
            {
                return ("Windows Operating System Version:  " + Environment.OSVersion.Version.ToString() + "\n");
            }
            catch
            {
                return "Windows Operating System Version could not be determined.\n";
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

                return ("Installed RAM:     " + (Capacity / 1073741824).ToString() +"\n");
            }
            catch
            {
                return "Amount of installed RAM could not be determined.\n";
            }
        }

        private static string ProcessorCount()
        {
            try
            {
                return ("ProcessorCount:  " + Environment.ProcessorCount.ToString() + "\n");
            }
            catch
            {
                return "The number of processors could not be determined.\n";
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
