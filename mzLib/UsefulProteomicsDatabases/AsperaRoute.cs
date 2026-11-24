using System;
using System.Threading;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The state of the Aspera route that every <see cref="PrideDownloadPolicy"/> in the process shares: how many
    /// sessions are running, when the last one started, and whether the route is set aside after a route failure.
    /// </summary>
    /// <remarks>
    /// One instance, <see cref="Shared"/>, serves the process, so two clients, or two policies, still stagger their
    /// starts and count against one cap: EBI sees one machine. Tests give a client its own instance and clock.
    /// </remarks>
    internal sealed class AsperaRoute
    {
        /// <summary>The process's Aspera route.</summary>
        internal static readonly AsperaRoute Shared = new(() => DateTimeOffset.UtcNow);

        private readonly object _lock = new();
        private readonly Func<DateTimeOffset> _clock;
        private int _running;
        private DateTimeOffset _lastStart = DateTimeOffset.MinValue;
        private DateTimeOffset _pausedUntil = DateTimeOffset.MinValue;
        private TaskCompletionSource _released = NewSignal();

        /// <param name="clock">Now, in UTC. Tests pass a clock they move by hand.</param>
        internal AsperaRoute(Func<DateTimeOffset> clock) => _clock = clock;

        /// <summary>Whether the route is set aside now.</summary>
        internal bool IsPaused
        {
            get
            {
                lock (_lock)
                    return _clock() < _pausedUntil;
            }
        }

        /// <summary>Sets the route aside for <paramref name="duration"/> from now, unless it already is for longer.</summary>
        internal void Pause(TimeSpan duration)
        {
            if (duration <= TimeSpan.Zero)
                return;
            lock (_lock)
            {
                DateTimeOffset until = _clock() + duration;
                if (until > _pausedUntil)
                    _pausedUntil = until;
            }
        }

        /// <summary>
        /// Waits until fewer than <paramref name="cap"/> sessions are running and <paramref name="startGap"/> has passed
        /// since the last start, then counts this session as started. Dispose the result when the session ends.
        /// </summary>
        /// <param name="cap">The most sessions at once.</param>
        /// <param name="startGap">The least time since the last start.</param>
        /// <param name="delay">Waits out the gap; the client's backoff delay, so tests need not sleep.</param>
        /// <param name="cancellationToken">Ends the wait.</param>
        internal async Task<IDisposable> EnterAsync(int cap, TimeSpan startGap,
            Func<TimeSpan, CancellationToken, Task> delay, CancellationToken cancellationToken)
        {
            while (true)
            {
                Task released;
                TimeSpan wait;
                lock (_lock)
                {
                    released = _released.Task;
                    wait = TimeSpan.Zero;
                    if (_running < cap)
                    {
                        DateTimeOffset now = _clock();
                        wait = _lastStart + startGap - now;
                        if (_lastStart == DateTimeOffset.MinValue || wait <= TimeSpan.Zero)
                        {
                            _running++;
                            _lastStart = now;
                            return new Session(this);
                        }
                    }
                }

                if (wait > TimeSpan.Zero)
                    await delay(wait, cancellationToken).ConfigureAwait(false);
                else
                    await released.WaitAsync(cancellationToken).ConfigureAwait(false);
            }
        }

        private void Leave()
        {
            TaskCompletionSource released;
            lock (_lock)
            {
                _running--;
                released = _released;
                _released = NewSignal();
            }
            released.TrySetResult();
        }

        private static TaskCompletionSource NewSignal() => new(TaskCreationOptions.RunContinuationsAsynchronously);

        private sealed class Session : IDisposable
        {
            private AsperaRoute _route;

            internal Session(AsperaRoute route) => _route = route;

            public void Dispose() => Interlocked.Exchange(ref _route, null)?.Leave();
        }
    }
}
