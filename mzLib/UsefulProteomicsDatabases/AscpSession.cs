using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One running <c>ascp</c>, as <see cref="PrideArchiveClient.DownloadFileViaAsperaAsync"/> watches it. Tests replace
    /// it (<see cref="PrideArchiveClient.StartAscp"/>) to script what ascp writes and says, with no ascp installed.
    /// </summary>
    internal interface IAscpSession : IDisposable
    {
        /// <summary>Completes with ascp's exit code once it has exited and its output has been read to the end.</summary>
        Task<int> Exited { get; }

        /// <summary>The last non-blank line ascp wrote, to either stream, or an empty string.</summary>
        string LastLine { get; }

        /// <summary>Stops ascp and every process it started. Does nothing once it has exited.</summary>
        void Kill();
    }

    /// <summary>The real <see cref="IAscpSession"/>: an ascp process, its output read as it arrives.</summary>
    internal sealed class AscpProcess : IAscpSession
    {
        private readonly Process _process;
        private string _lastLine = string.Empty;

        private AscpProcess(Process process)
        {
            _process = process;
            _process.OutputDataReceived += OnLine;
            _process.ErrorDataReceived += OnLine;
        }

        /// <summary>Starts ascp.</summary>
        /// <exception cref="Win32Exception">The executable could not be started (missing, not executable).</exception>
        /// <exception cref="InvalidOperationException">The process did not start.</exception>
        internal static IAscpSession Start(ProcessStartInfo startInfo)
        {
            var process = new Process { StartInfo = startInfo };
            var session = new AscpProcess(process);
            try
            {
                if (!process.Start())
                    throw new InvalidOperationException("The process did not start.");
                process.BeginOutputReadLine();
                process.BeginErrorReadLine();
            }
            catch
            {
                process.Dispose();
                throw;
            }

            session.Exited = session.WaitForExitAsync();
            return session;
        }

        public Task<int> Exited { get; private set; }

        public string LastLine => Volatile.Read(ref _lastLine);

        public void Kill()
        {
            try
            {
                _process.Kill(entireProcessTree: true);
            }
            catch (InvalidOperationException)
            {
                // Already exited.
            }
            catch (Win32Exception)
            {
                // Exiting, or a child the OS no longer lets us end; WaitForExitAsync still reports the exit.
            }
        }

        public void Dispose() => _process.Dispose();

        private async Task<int> WaitForExitAsync()
        {
            // Waits for the redirected streams to be read to the end too, so LastLine is final when this completes.
            await _process.WaitForExitAsync(CancellationToken.None).ConfigureAwait(false);
            return _process.ExitCode;
        }

        private void OnLine(object sender, DataReceivedEventArgs e)
        {
            if (!string.IsNullOrWhiteSpace(e.Data))
                Volatile.Write(ref _lastLine, e.Data.Trim());
        }
    }
}
