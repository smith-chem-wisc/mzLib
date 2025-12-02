using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;

namespace IsobaricQuant
{
    public class IsobaricQuantEngine
    {
        private List<(int peptideFullSequenceHash, int[] reporterIntensities)> theInput;
        

        public IsobaricQuantEngine(List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> myFractionInput)
        {
            theInput = CollectPsmsFromFractions(myFractionInput);
        }
        public ConcurrentDictionary<int, ConcurrentDictionary<int, int>> Process(AggregateType aggregateType, int topN = int.MaxValue, double lowFraction = 0.05, bool normalize = true, int? referenceChannel = null)
        {
            var k = GroupAllRawPsmChannelIntensities(theInput);
            var l = CombineGroupedPsmIntensitiesToSingleValue(k, aggregateType, topN);
            var m = RemoveLowIntensityPsms(l, lowFraction);
            if (normalize)
            {
                if (referenceChannel == null)
                {
                    var n = DetermineChannelMedians(m);
                    var o = ComputeChannelShifts(n);
                    return ApplyChannelShifts(m, o);
                }

                return SubtractReferenceChannel(m, referenceChannel);
            }
            else return m;
        }

        //public IsobaricQuantResults Run()
        //{
        //    Process();
        //}

        /// <summary>
        /// This could be the input of values from PSMs of a single file. On the input side, the peptideFullSequenceHash is exactly that and it may be 
        /// duplicated in the list of psmReporterValues. This is psm level data. So, there should be one copy of each peptideFullSequenceHash for each
        /// PSM with the same full sequence. The int[] reporterIntensities will be all reporter ion intensities for the corresponding psm. These can be made ints by
        /// rounding or multiplying by a large number and rounding, etc. This is done for speed and efficiency
        /// What you get in return is an aggregated list of reporter intensities for each full sequence and each channel
        /// </summary>
        /// <param name="psmReporterValues"></param>
        /// <returns></returns>
        private ConcurrentDictionary<int, ConcurrentDictionary<int, ConcurrentBag<int>>> GroupAllRawPsmChannelIntensities(List<(int peptideFullSequenceHash, int[] reporterIntensities)> psmReporterValues)
        {
            var myOut = new ConcurrentDictionary<int, ConcurrentDictionary<int, ConcurrentBag<int>>>();

            if (psmReporterValues == null || psmReporterValues.Count == 0)
            {
                return myOut;
            }

            int expectedChannels = psmReporterValues[0].reporterIntensities.Length;

            // Parallel aggregation
            Parallel.ForEach(psmReporterValues, tuple =>
            {
                // Deconstruct for clarity
                var (peptideFullSequenceHash, intensities) = tuple;

                // Skip malformed entries
                if (intensities == null || intensities.Length != expectedChannels)
                    return;

                // Get or create inner dictionary for this plex id
                var channelDict = myOut.GetOrAdd(peptideFullSequenceHash, _ => new ConcurrentDictionary<int, ConcurrentBag<int>>());

                // Add each channel's intensity
                for (int channelIndex = 0; channelIndex < intensities.Length; channelIndex++)
                {
                    var bag = channelDict.GetOrAdd(channelIndex, _ => new ConcurrentBag<int>());
                    bag.Add(intensities[channelIndex]);
                }
            });

            return myOut;
        }
        /// <summary>
        /// Aggregates peptide intensity values by peptide and channel using the specified aggregation method.
        /// </summary>
        /// <remarks>For each peptide and channel, the aggregation is performed as follows: Max returns
        /// the maximum intensity, Median returns the median intensity, and SumTopN returns the sum of the top N
        /// intensities. If a channel has no intensity values, the aggregated value is 0. The method is thread-safe and
        /// optimized for concurrent execution.</remarks>
        /// <param name="peptideIntensitiesByPlexAndFile">A concurrent dictionary mapping peptide identifiers to channel-specific collections of intensity values.
        /// Each inner dictionary maps a channel index to a collection of intensity values for that peptide and channel.</param>
        /// <param name="aggregateType">The aggregation method to apply to the intensity values for each peptide and channel. Supported values
        /// include Max, Median, and SumTopN.</param>
        /// <param name="topN">The maximum number of top intensity values to include when using the SumTopN aggregation method. Ignored for
        /// other aggregation types. Must be greater than 0 to have an effect.</param>
        /// <returns>A concurrent dictionary mapping each peptide identifier to a dictionary of aggregated intensity values by
        /// channel. If the input is null or empty, returns an empty dictionary.</returns>
        private ConcurrentDictionary<int, ConcurrentDictionary<int, int>> CombineGroupedPsmIntensitiesToSingleValue(ConcurrentDictionary<int, ConcurrentDictionary<int, ConcurrentBag<int>>> peptideIntensitiesByPlexAndFile, AggregateType aggregateType, int topN = int.MaxValue)
        {
            var myOut = new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();
            if (peptideIntensitiesByPlexAndFile == null || peptideIntensitiesByPlexAndFile.IsEmpty)
                return myOut;

            Parallel.ForEach(peptideIntensitiesByPlexAndFile, outerKvp =>
            {
                var peptideHash = outerKvp.Key;
                var channelBags = outerKvp.Value;
                var aggregatedChannels = new ConcurrentDictionary<int, int>();

                foreach (var channelKvp in channelBags)
                {
                    int channelIndex = channelKvp.Key;
                    var bag = channelKvp.Value;

                    // Empty bag -> 0
                    if (bag == null)
                    {
                        aggregatedChannels[channelIndex] = 0;
                        continue;
                    }

                    int[] values = bag.ToArray();
                    if (values.Length == 0)
                    {
                        aggregatedChannels[channelIndex] = 0;
                        continue;
                    }

                    int result;
                    switch (aggregateType)
                    {
                        case AggregateType.Max:
                            {
                                int max = int.MinValue;
                                for (int i = 0; i < values.Length; i++)
                                {
                                    if (values[i] > max) max = values[i];
                                }
                                result = max;
                                break;
                            }
                        case AggregateType.Median:
                            {
                                Array.Sort(values); // ascending
                                int len = values.Length;
                                if ((len & 1) == 1)
                                {
                                    result = values[len / 2];
                                }
                                else
                                {
                                    // Truncated average of middle two (deterministic int)
                                    int a = values[(len / 2) - 1];
                                    int b = values[len / 2];
                                    result = (a + b) / 2;
                                }
                                break;
                            }
                        case AggregateType.SumTopN:
                            {
                                if (topN <= 0)
                                {
                                    result = 0;
                                    break;
                                }

                                if (values.Length <= topN)
                                {
                                    long sumAll = 0;
                                    for (int i = 0; i < values.Length; i++)
                                        sumAll += values[i];
                                    // Clamp to int if overflow (unlikely for reporter ints, but safe)
                                    result = sumAll > int.MaxValue ? int.MaxValue : (int)sumAll;
                                    break;
                                }

                                // Maintain a min-heap of size topN for best performance when values.Length >> topN
                                // Heap stored in-place in an array 'heap'
                                int[] heap = new int[topN];
                                int size = 0;

                                void HeapifyDown(int idx)
                                {
                                    while (true)
                                    {
                                        int left = (idx * 2) + 1;
                                        if (left >= size) break;
                                        int right = left + 1;
                                        int smallest = left;
                                        if (right < size && heap[right] < heap[left]) smallest = right;
                                        if (heap[smallest] >= heap[idx]) break;
                                        (heap[idx], heap[smallest]) = (heap[smallest], heap[idx]);
                                        idx = smallest;
                                    }
                                }

                                // Build heap from initial topN values
                                for (int i = 0; i < values.Length; i++)
                                {
                                    int v = values[i];
                                    if (size < topN)
                                    {
                                        heap[size++] = v;
                                        if (size == topN)
                                        {
                                            for (int h = (size / 2) - 1; h >= 0; h--)
                                                HeapifyDown(h);
                                        }
                                    }
                                    else if (v > heap[0])
                                    {
                                        heap[0] = v;
                                        HeapifyDown(0);
                                    }
                                }

                                long sum = 0;
                                for (int i = 0; i < size; i++)
                                    sum += heap[i];

                                result = sum > int.MaxValue ? int.MaxValue : (int)sum;
                                break;
                            }
                        default:
                            result = 0;
                            break;
                    }

                    aggregatedChannels[channelIndex] = result;
                }

                myOut[peptideHash] = aggregatedChannels;
            });
            return myOut;
        }
        /// <summary>
        /// Removes the bottom fraction of peptides ranked by total intensity (sum across channels).
        /// </summary>
        /// <param name="peptideIntensitiesAggregated">
        /// Outer key = peptide id; inner dictionary maps channel index → aggregated intensity for that peptide.
        /// </param>
        /// <param name="fraction">
        /// Fraction in [0,1]. Removes floor(fraction * count) peptides with the smallest total intensity.
        /// ≤ 0 removes nothing; ≥ 1 removes all.
        /// </param>
        /// <returns>
        /// A new ConcurrentDictionary containing only the retained peptides. Returns empty for null/empty input.
        /// </returns>
        /// <remarks>
        /// Implementation:
        /// - Snapshots the input for stable enumeration.
        /// - Computes per-peptide sums in parallel (long accumulator to avoid overflow).
        /// - Selects the k smallest sums (k = floor(fraction * N)) using an O(N log k) max-heap (avoids full sort).
        /// - Builds and returns a new dictionary excluding those k peptides. Input is not mutated (thread-safe).
        /// Complexity: O(P + N log k) time, O(N) memory; P = total channel entries enumerated, N = peptide count.
        /// Ties at the cutoff are resolved by heap ordering (not stable across differing input orders).
        /// </remarks>
        private ConcurrentDictionary<int, ConcurrentDictionary<int, int>> RemoveLowIntensityPsms(ConcurrentDictionary<int, ConcurrentDictionary<int, int>> peptideIntensitiesAggregated, double fraction = 0.05)
        {
            // Fast exits
            if (peptideIntensitiesAggregated == null || peptideIntensitiesAggregated.IsEmpty)
                return new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();

            int total = peptideIntensitiesAggregated.Count;
            if (fraction <= 0) return peptideIntensitiesAggregated; // nothing to remove
            if (fraction >= 1) return new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>(); // remove all

            int removeCount = (int)Math.Floor(fraction * total);
            if (removeCount <= 0) return peptideIntensitiesAggregated;
            if (removeCount >= total) return new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();

            // Snapshot for consistent, indexable enumeration
            var snapshot = peptideIntensitiesAggregated.ToArray(); // KeyValuePair<int, ConcurrentDictionary<int,int>>[]
            var sums = new (int key, long sum)[snapshot.Length];

            // Parallel sum per peptide (use long to avoid overflow)
            Parallel.For(0, snapshot.Length, i =>
            {
                var kvp = snapshot[i];
                long s = 0;
                foreach (var v in kvp.Value.Values)
                    s += v;
                sums[i] = (kvp.Key, s);
            });

            // Select the removeCount smallest sums using a fixed-size max-heap (O(N log k))
            int k = removeCount;
            var heap = new (long sum, int key)[k];
            int size = 0;

            void HeapifyDown(int idx)
            {
                while (true)
                {
                    int left = (idx * 2) + 1;
                    if (left >= size) break;
                    int right = left + 1;
                    int largest = left;
                    if (right < size && heap[right].sum > heap[left].sum) largest = right;
                    if (heap[largest].sum <= heap[idx].sum) break;
                    (heap[idx], heap[largest]) = (heap[largest], heap[idx]);
                    idx = largest;
                }
            }

            // Build/maintain max-heap of the k smallest sums
            for (int i = 0; i < sums.Length; i++)
            {
                var (key, sum) = sums[i];
                if (size < k)
                {
                    heap[size++] = (sum, key);
                    if (size == k)
                    {
                        for (int h = (size / 2) - 1; h >= 0; h--)
                            HeapifyDown(h);
                    }
                }
                else if (sum < heap[0].sum)
                {
                    heap[0] = (sum, key);
                    HeapifyDown(0);
                }
            }

            // Keys to remove (bottom fraction by total intensity)
            var removeKeys = new HashSet<int>();
            for (int i = 0; i < size; i++)
                removeKeys.Add(heap[i].key);

            // Build result with retained entries
            var myOut = new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();
            for (int i = 0; i < snapshot.Length; i++)
            {
                var kvp = snapshot[i];
                if (!removeKeys.Contains(kvp.Key))
                    myOut.TryAdd(kvp.Key, kvp.Value);
            }
            return myOut;
        }
        private int[] DetermineChannelMedians(ConcurrentDictionary<int, ConcurrentDictionary<int, int>> rawValues)
        {
            // Handle null/empty input
            if (rawValues == null || rawValues.IsEmpty)
                return Array.Empty<int>();

            // Snapshot and find max channel key
            var snapshot = rawValues.ToArray();
            int maxKey = -1;
            for (int i = 0; i < snapshot.Length; i++)
            {
                var channels = snapshot[i].Value;
                foreach (var kvp in channels)
                {
                    if (kvp.Key > maxKey) maxKey = kvp.Key;
                }
            }

            if (maxKey < 0)
                return Array.Empty<int>();

            // Length is max key + 1 so that index == channel key
            int length = maxKey + 1;
            var perChannelValues = new List<int>[length];
            var locks = new object[length];
            for (int i = 0; i < length; i++)
            {
                perChannelValues[i] = new List<int>(256); // small default capacity
                locks[i] = new object();
            }

            // Aggregate values per channel (lock-stripe per index; keys < 64 keeps contention low)
            for (int i = 0; i < snapshot.Length; i++)
            {
                var channels = snapshot[i].Value;
                foreach (var kvp in channels)
                {
                    int channel = kvp.Key;
                    int value = kvp.Value;
                    lock (locks[channel])
                    {
                        perChannelValues[channel].Add(value);
                    }
                }
            }

            // Compute medians; parallel over channels for speed
            var medians = new int[length];
            Parallel.For(0, length, ch =>
            {
                var values = perChannelValues[ch];
                if (values.Count == 0)
                {
                    medians[ch] = 0;
                    return;
                }

                values.Sort(); // ascending
                int count = values.Count;
                if ((count & 1) == 1)
                {
                    medians[ch] = values[count / 2];
                }
                else
                {
                    int a = values[(count / 2) - 1];
                    int b = values[count / 2];
                    medians[ch] = (a + b) / 2; // integer median (truncated average of middle two)
                }
            });

            return medians;
        }
        /// <summary>
        /// Computes per-channel shifts by centering the supplied channel medians on their overall median.
        /// </summary>
        /// <param name="channelMedians">
        /// Array of per-channel median intensities. Must be non-null; length typically < 64.
        /// </param>
        /// <returns>
        /// An int[] of the same length where each element is (channelMedians[i] - overallMedian).
        /// Returns an empty array if input is null or length == 0.
        /// </returns>
        /// <remarks>
        /// Overall median is computed on a sorted copy:
        /// - Odd length: middle element.
        /// - Even length: truncated average of the two middle elements.
        /// The operation preserves integer arithmetic (no rounding issues).
        /// Time: O(C log C) due to sort (C = channel count, small); Space: O(C) for the copy.
        /// </remarks>
        private int[] ComputeChannelShifts(int[] channelMedians)
        {
            // Null or empty -> nothing to shift
            if (channelMedians == null || channelMedians.Length == 0)
                return Array.Empty<int>();

            // Copy and sort to compute overall median efficiently (channel count < 64)
            int[] tmp = new int[channelMedians.Length];
            Array.Copy(channelMedians, tmp, tmp.Length);
            Array.Sort(tmp);

            int len = tmp.Length;
            int overallMedian;
            if ((len & 1) == 1)
            {
                overallMedian = tmp[len / 2];
            }
            else
            {
                int a = tmp[(len / 2) - 1];
                int b = tmp[len / 2];
                overallMedian = (a + b) / 2; // truncated integer midpoint
            }

            // Produce shifts: value - overallMedian
            int[] shifts = new int[channelMedians.Length];
            for (int i = 0; i < channelMedians.Length; i++)
            {
                shifts[i] = channelMedians[i] - overallMedian;
            }

            return shifts;
        }
        /// <summary>
        /// Applies per-channel integer shifts to all peptide intensities in-place.
        /// For each outer entry (peptide), every inner (channel,value) pair where channel < channelShifts.Length
        /// is updated to (value + channelShifts[channel]). Returns the (now shifted) input dictionary.
        /// </summary>
        /// <param name="rawIntensities">Outer: peptide id; Inner: channel index → intensity.</param>
        /// <param name="channelShifts">Per-channel shift values; index corresponds to channel key.</param>
        /// <returns>The same ConcurrentDictionary instance with updated values (no new allocations beyond loop work). If shifts are null/empty, input is returned unchanged.</returns>
        private ConcurrentDictionary<int, ConcurrentDictionary<int, int>> ApplyChannelShifts(
            ConcurrentDictionary<int, ConcurrentDictionary<int, int>> rawIntensities,
            int[] channelShifts)
        {
            if (rawIntensities == null || rawIntensities.IsEmpty)
                return rawIntensities ?? new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();
            if (channelShifts == null || channelShifts.Length == 0)
                return rawIntensities;

            // Parallel over peptides; inner enumeration is a snapshot so concurrent writes are safe.
            Parallel.ForEach(rawIntensities, peptideEntry =>
            {
                var inner = peptideEntry.Value;
                foreach (var channelEntry in inner)
                {
                    int channel = channelEntry.Key;
                    if ((uint)channel >= (uint)channelShifts.Length) // unsigned compare avoids negative issues
                        continue;

                    int shift = channelShifts[channel];
                    if (shift == 0)
                        continue;

                    int original = channelEntry.Value;
                    long newValLong = (long)original + shift; // promote to avoid overflow
                    int newVal = newValLong > int.MaxValue
                        ? int.MaxValue
                        : newValLong < int.MinValue
                            ? int.MinValue
                            : (int)newValLong;

                    // Update in-place (indexer is thread-safe for ConcurrentDictionary)
                    inner[channel] = newVal;
                }
            });

            return rawIntensities;
        }

        private ConcurrentDictionary<int, ConcurrentDictionary<int, int>> SubtractReferenceChannel(ConcurrentDictionary<int, ConcurrentDictionary<int, int>> m, int? referenceChannel)
        {
            // Validate inputs
            if (m == null || m.IsEmpty)
                return m ?? new ConcurrentDictionary<int, ConcurrentDictionary<int, int>>();
            if (referenceChannel is null || referenceChannel < 0)
                return m;

            int refCh = referenceChannel.Value;

            // Parallelize across peptides; subtract the per-peptide reference-channel value from all its channels
            Parallel.ForEach(m, peptideEntry =>
            {
                var channels = peptideEntry.Value;

                // Try get the reference value for this peptide; if missing, skip subtraction for this peptide
                if (!channels.TryGetValue(refCh, out int refValue))
                    return;

                if (refValue == 0)
                    return; // nothing to change

                // Subtract with overflow protection
                foreach (var kv in channels)
                {
                    int channel = kv.Key;
                    int val = kv.Value;

                    long diff = (long)val - refValue;
                    int newVal = diff > int.MaxValue
                        ? int.MaxValue
                        : diff < int.MinValue
                            ? int.MinValue
                            : (int)diff;

                    channels[channel] = newVal;
                }
            });

            return m;
        }
        private List<(int peptideFullSequenceHash, int[] reporterIntensities)> CollectPsmsFromFractions(List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> myFractionInput)
        {
            // Fast exits
            if (myFractionInput == null || myFractionInput.Count == 0)
            {
                new List<(int peptideFullSequenceHash, int[] reporterIntensities)>();
            }

            // Determine expected channel count from the first valid entry (for consistency)
            int? expectedChannels = null;
            for (int i = 0; i < myFractionInput.Count && expectedChannels is null; i++)
            {
                var frac = myFractionInput[i];
                if (frac == null) continue;
                for (int j = 0; j < frac.Count; j++)
                {
                    var arr = frac[j].reporterIntensities;
                    if (arr != null)
                    {
                        expectedChannels = arr.Length;
                        break;
                    }
                }
            }

            // If no valid arrays were found, result is empty
            if (expectedChannels is null)
            {
                return new List<(int peptideFullSequenceHash, int[] reporterIntensities)>(0);
            }

            // Collect in a thread-safe sink, then materialize to List at the end
            var sink = new ConcurrentBag<(int peptideFullSequenceHash, int[] reporterIntensities)>();

            // Parallelize across fractions for maximum throughput
            Parallel.ForEach(myFractionInput, fraction =>
            {
                if (fraction == null || fraction.Count == 0)
                    return;

                // Iterate fraction entries; filter fast
                for (int idx = 0; idx < fraction.Count; idx++)
                {
                    var entry = fraction[idx];
                    var arr = entry.reporterIntensities;
                    if (arr == null || arr.Length != expectedChannels.Value)
                        continue; // skip malformed/inconsistent entries

                    // Add to sink (ConcurrentBag is lock-free for adds)
                    sink.Add(entry);
                }
            });

            // Materialize result; pre-size list to avoid reallocations
            theInput = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>(sink.Count);
            foreach (var tuple in sink)
            {
                theInput.Add(tuple);
            }
            return theInput;
        }
        public enum AggregateType
        {
            SumTopN,
            Max,
            Median
        }
    }
}
