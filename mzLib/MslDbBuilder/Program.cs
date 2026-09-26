using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Threading;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Chemistry;
using Chromatography.RetentionTimePrediction;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using UsefulProteomicsDatabases;

// -------------------------------------------------------------------------------------------------
// MslDbBuilder — converts protein FASTA databases into MetaMorpheus .msl spectral libraries that
// hold PREFABRICATED THEORETICAL FRAGMENTATION (a fragmentation cache for the ManySearch GPU path).
//
// For each database: load targets + decoys (Reverse), digest and fragment exactly as the search
// would, and write one .msl per database. Each library entry stores a peptide's theoretical
// fragment-ion masses (intensities are placeholders — MetaMorpheus scoring uses EXPERIMENTAL
// intensities, so the library's intensities are intentionally unused).
//
// Digestion/mods: trypsin, 0 MISSED CLEAVAGES (the library defines the search peptide space),
// length >=7, variable initiator-Met, fixed Carbamidomethyl on C, variable Oxidation on M,
// max 2 mods. Decoys: DecoyType.Reverse.
//
// Usage:  MslDbBuilder <outputDir> <fasta1> [fasta2 ...]
// -------------------------------------------------------------------------------------------------

// Curated "Common Fixed"/"Common Variable" list — has "Carbamidomethyl on C", "Oxidation on M"
// with ID already == IdWithMotif (matches the toml's ListOfMods names exactly).
const string ModsFile = @"E:\CodeReview\metaProteomics\MetaMorpheus\MetaMorpheus\EngineLayer\Mods\aListOfmods.txt";
const int NominalCharge = 2;                  // precursor charge stored per entry (fragment masses are charge-independent)
const DissociationType Dissociation = DissociationType.HCD;

bool verify = args.Contains("--verify");
bool readMode = args.Contains("--read");
// --lean: build a fragment-LESS index (precursor mass + charge + iRT + sequence + accession, NO stored
// fragments). The search fragments candidates on the fly (cheap), so storing fragments only bloats the
// .msl ~10x and makes the search disk-bound on fragment fetches. Lean files are ~10x smaller and faster
// to both build and read.
bool lean = args.Contains("--lean");
var positional = args.Where(a => a != "--verify" && a != "--read" && a != "--rtcalib" && a != "--checkrt" && a != "--lean" && a != "--merge" && a != "--shardbuild" && a != "--streamtest" && a != "--probe" && a != "--noirt").ToArray();

// --streamtest <N> [outPath]: isolate MslWriter.WriteStreaming. Generates N synthetic LEAN entries
// (empty fragment list, unique sequence + "db|acc" per entry => big string table), writes via
// WriteStreaming, reads back with LoadIndexOnly, and reports header-offset consistency + entry count.
if (args.Contains("--streamtest"))
{
    int n = positional.Length >= 1 ? int.Parse(positional[0]) : 1_000_000;
    string outPath = positional.Length >= 2 ? positional[1] : @"E:\CodeReview\metaProteomics\streamtest.msl";
    Console.WriteLine($"--streamtest: writing {n} lean entries -> {outPath}");

    IEnumerable<MslLibraryEntry> Gen()
    {
        var aa = "ACDEFGHIKLMNPQRSTVWY";
        for (int i = 0; i < n; i++)
        {
            // Build a deterministic, unique-ish sequence so the string table grows like the real shard.
            int len = 7 + (i % 20);
            var sb = new System.Text.StringBuilder(len);
            int x = i;
            for (int k = 0; k < len; k++) { sb.Append(aa[x % 20]); x = x * 1103515245 + 12345; x &= 0x7fffffff; }
            string seq = sb.ToString();
            var spec = new LibrarySpectrum(seq, 400.0 + (i % 100000) * 0.01, 2,
                new List<MatchedFragmentIon>(), rt: (i % 90) + 0.5, isDecoy: (i & 1) == 1);
            var e = MslLibraryEntry.FromLibrarySpectrum(spec);
            e.ProteinAccession = "db" + (i % 50) + "|ACC" + i;
            yield return e;
        }
    }

    if (n >= 0)
        MslWriter.WriteStreaming(outPath, Gen());
    long fsz = new FileInfo(outPath).Length;
    Console.WriteLine($"  wrote {fsz} bytes");

    using (var lib = MslLibrary.LoadIndexOnly(outPath))
    {
        int qcount = 0;
        foreach (var _ in lib.QueryMzWindow(float.MinValue, float.MaxValue)) qcount++;
        Console.WriteLine($"  LoadIndexOnly OK. QueryMzWindow(all) count = {qcount} (expected {n})");
        if (qcount > 0)
        {
            var e0 = lib.GetEntry(0);
            Console.WriteLine($"  GetEntry(0): seq='{e0?.FullSequence}' acc='{e0?.ProteinAccession}' mz={e0?.PrecursorMz}");
        }
    }
    return 0;
}

// --probe <path>: open an EXISTING .msl with both LoadIndexOnly and Load, report counts + full exceptions.
if (args.Contains("--probe"))
{
    string p = positional[0];
    Console.WriteLine($"--probe: {p}  size={new FileInfo(p).Length}");
    try
    {
        using var lib = MslLibrary.LoadIndexOnly(p);
        int q = 0; float mn = float.MaxValue, mx = float.MinValue;
        foreach (var e in lib.QueryMzWindow(float.MinValue, float.MaxValue)) { q++; if (e.PrecursorMz < mn) mn = e.PrecursorMz; if (e.PrecursorMz > mx) mx = e.PrecursorMz; }
        Console.WriteLine($"  LoadIndexOnly: QueryMzWindow(all) count={q} mzMin={mn} mzMax={mx}");
        int q2 = 0; foreach (var e in lib.QueryMzWindow(0f, 5000f)) q2++;
        Console.WriteLine($"  QueryMzWindow(0,5000) count={q2}");
    }
    catch (Exception ex) { Console.WriteLine($"  LoadIndexOnly THREW: {ex.GetType().Name}: {ex.Message}\n{ex.StackTrace}"); }
    try
    {
        using var lib = MslLibrary.Load(p);
        int c = lib.GetAllEntries(true).Count();
        Console.WriteLine($"  Load: GetAllEntries(true) count={c}");
    }
    catch (Exception ex) { Console.WriteLine($"  Load THREW: {ex.GetType().Name}: {ex.Message}"); }
    return 0;
}

// --- load the exact MetaMorpheus mods (needed by both modes) ---
var allMods = PtmListLoader.ReadModsFromFile(ModsFile, out _).ToList();

// --read mode: reconstruct searchable peptides from .msl and confirm double-precision
// refragmentation matches the stored fragmentation (the foundation for the library-based search).
if (readMode)
{
    if (positional.Length < 1)
    {
        Console.Error.WriteLine("Usage: MslDbBuilder --read <lib1.msl> [lib2.msl ...]");
        return 1;
    }
    var allKnownMods = new Dictionary<string, Modification>();
    foreach (var m in allMods) allKnownMods[m.IdWithMotif] = m; // last wins
    foreach (var msl in positional)
        ReadAndCheckPeptides(msl, allKnownMods);
    return 0;
}

// --rtcalib mode: measure how accurately an RT predictor places confident base-search PSMs on THIS
// gradient. Regress observed RT vs predicted iRT and report the residual SD — that residual is the RT
// window width, hence the selectivity multiplier for an RT-indexed .msl candidate filter.
if (args.Contains("--rtcalib"))
{
    if (positional.Length < 1)
    {
        Console.Error.WriteLine("Usage: MslDbBuilder --rtcalib <BasePSMs.psmtsv>");
        return 1;
    }
    var allKnownMods = new Dictionary<string, Modification>();
    foreach (var m in allMods) allKnownMods[m.IdWithMotif] = m;
    // Default Chronologer (accurate); pass "ssrcalc" as 2nd positional for the fast analytical predictor.
    var predType = positional.Length > 1 && positional[1].Equals("ssrcalc", StringComparison.OrdinalIgnoreCase)
        ? PredictorType.SSRCalc3 : PredictorType.Chronologer;
    RunRtCalib(positional[0], allKnownMods, predType);
    return 0;
}

// --checkrt mode: confirm stored iRT landed in the index (LoadIndexOnly reads precursor metadata incl. Irt).
if (args.Contains("--checkrt"))
{
    foreach (var msl in positional)
    {
        using var lib = MslLibrary.LoadIndexOnly(msl);
        var irts = new List<float>();
        foreach (var e in lib.QueryMzWindow(float.MinValue, float.MaxValue)) irts.Add(e.Irt);
        int nonzero = irts.Count(x => x != 0f);
        Console.WriteLine($"CHECKRT {Path.GetFileName(msl)}: entries={irts.Count} nonZeroIrt={nonzero} " +
            $"min={(irts.Count > 0 ? irts.Min() : 0):F2} max={(irts.Count > 0 ? irts.Max() : 0):F2}");
    }
    return 0;
}

// --merge mode: combine many (lean) .msl into ONE merged .msl, stamping each entry's accession with its
// source proteome ("<db>|<accession>") so per-database results can still be recovered. A single merged
// file = a single LoadIndexOnly open instead of thousands of random file opens (the search bottleneck).
//   MslDbBuilder --merge <outMsl> <in1.msl> [in2.msl ...]   (or  --merge <outMsl> @<listfile>)
if (args.Contains("--merge"))
{
    if (positional.Length < 2)
    {
        Console.Error.WriteLine("Usage: MslDbBuilder --merge <outMsl> <in1.msl ...|@listfile>");
        return 1;
    }
    string mergedOut = positional[0];
    string[] inputs = positional.Length == 2 && positional[1].StartsWith("@")
        ? File.ReadAllLines(positional[1].Substring(1)).Select(l => l.Trim()).Where(l => l.Length > 0).ToArray()
        : positional.Skip(1).ToArray();

    var merged = new List<MslLibraryEntry>(1 << 20);
    int done = 0;
    foreach (var inPath in inputs)
    {
        string db = Path.GetFileNameWithoutExtension(inPath);
        using (var lib = MslLibrary.Load(inPath))
        {
            foreach (var e in lib.GetAllEntries(includeDecoys: true))
            {
                e.ProteinAccession = db + "|" + (e.ProteinAccession ?? "");
                merged.Add(e);
            }
        }
        if (++done % 1000 == 0) Console.WriteLine($"  merged {done}/{inputs.Length} ({merged.Count} entries)");
    }
    Console.WriteLine($"Writing merged .msl: {merged.Count} entries from {inputs.Length} databases -> {mergedOut}");
    MslLibrary.Save(mergedOut, merged);
    Console.WriteLine($"Done. {new FileInfo(mergedOut).Length / 1048576} MB");
    return 0;
}

// --shardbuild mode: build N LEAN merged shards DIRECTLY from FASTAs (no per-db files), in PARALLEL,
// streaming so entries are never all held in RAM. Each shard is one MslWriter.WriteStreaming over an
// IEnumerable that loads/digests/iRT-predicts one database at a time and yields its (accession-stamped
// "db|acc") lean entries. Databases are kept WHOLE within a shard (so the search can union per-db groups
// across shards) and balanced by FASTA size. A single .msl maxes at int.MaxValue (~2.15B) entries, so the
// full bacteria set (~3.46B) MUST be sharded.
//   MslDbBuilder --shardbuild <outDir> <numShards> <fasta...|@listfile>
if (args.Contains("--shardbuild"))
{
    if (positional.Length < 3)
    {
        Console.Error.WriteLine("Usage: MslDbBuilder --shardbuild <outDir> <numShards> <fasta...|@listfile>");
        return 1;
    }
    string shardOutDir = positional[0];
    Directory.CreateDirectory(shardOutDir);
    int numShards = int.Parse(positional[1]);
    string[] allFastas = positional.Length == 3 && positional[2].StartsWith("@")
        ? File.ReadAllLines(positional[2].Substring(1)).Select(l => l.Trim()).Where(l => l.Length > 0).ToArray()
        : positional.Skip(2).ToArray();

    // Balance databases across shards by FASTA size (greedy: biggest first -> least-loaded shard),
    // keeping each database whole in one shard.
    var shardFastas = new List<string>[numShards];
    for (int i = 0; i < numShards; i++) shardFastas[i] = new List<string>();
    var shardBytes = new long[numShards];
    foreach (var f in allFastas.OrderByDescending(p => new FileInfo(p).Length))
    {
        int si = 0; for (int j = 1; j < numShards; j++) if (shardBytes[j] < shardBytes[si]) si = j;
        shardFastas[si].Add(f); shardBytes[si] += new FileInfo(f).Length;
    }

    var sFixedMods = new List<Modification> { GetMod("Carbamidomethyl on C") };
    var sVariableMods = new List<Modification> { GetMod("Oxidation on M") };
    var sDigest = new DigestionParams(maxMissedCleavages: 0);

    // SPEED/CRASH design: construct ONE Chronologer predictor (concurrent per-shard construction races
    // on the shared TorchSharp weights temp file and crashes). Chronologer inference already saturates
    // all cores internally, so iRT is the serial floor regardless of shard count — N predictors only
    // oversubscribe. We therefore serialize iRT inference behind a lock (free, since one predictor uses
    // every core) while digestion + shard WriteStreaming run in parallel across shards, streaming so
    // memory stays bounded (only the current db's peptides per worker).
    bool noIrt = args.Contains("--noirt");   // skip Chronologer entirely (rt:null) for a fast structural build
    int coreBudget = Math.Max(1, Environment.ProcessorCount - 2);
    // Cap how many shards build CONCURRENTLY, independent of shard COUNT. Fewer-concurrent + more-shards =>
    // shards finalize in WAVES: each completed .msl is a durable checkpoint, so an interruption loses only the
    // in-flight wave's digest+write (iRT is already banked in irt_cache.tsv), not the whole run. Also bounds
    // peak RAM (one db's peptides resident per concurrent shard). Override via MSLDB_SHARD_CONCURRENCY.
    int shardConcurrency = int.TryParse(Environment.GetEnvironmentVariable("MSLDB_SHARD_CONCURRENCY"), out var scov) && scov > 0
        ? Math.Min(numShards, scov) : Math.Min(numShards, 6);
    var irtLock = new object();
    long[] irtBusyTicks = new long[1];   // total wall the (serialized) Chronologer predictor spent actually inferring
    IRetentionTimePredictor sharedPred = noIrt ? null : RetentionTimePredictorFactory.Create(PredictorType.Chronologer);

    // GLOBAL iRT cache (Option A): iRT is a pure function of the peptide sequence, so a peptide that recurs
    // across organisms is predicted ONCE and reused everywhere. Hash-keyed (ConcurrentDictionary) => O(1)
    // membership/lookup regardless of size; no sorting, no binary search. Populated lazily as databases are
    // streamed (no upfront pass). null value = a sequence the predictor couldn't score (cached so it isn't
    // re-attempted). Persisted to irt_cache.tsv so a restart re-pays cheap digest+write but NEVER iRT.
    var irtCache = new ShardedSeqCache();
    string cachePath = Path.Combine(shardOutDir, "irt_cache.tsv");
    if (!noIrt && File.Exists(cachePath))
    {
        // Parallel preload — safe now that the cache is SHARDED: each line routes to an independent per-shard
        // ConcurrentDictionary, so concurrent inserts never race on one dict's resize (the bug that crashed an
        // earlier parallel attempt against the single dict). Single-threaded insert of ~900M entries took
        // >60 min with GC thrash; fanning the parse+insert across cores cuts it to a few minutes.
        Parallel.ForEach(File.ReadLines(cachePath), line =>
        {
            int t = line.IndexOf('\t');
            if (t < 0) return; // skip a partial trailing line from a prior crash
            string val = line.Substring(t + 1);
            irtCache[line.Substring(0, t)] = val.Length == 0 ? (double?)null : double.Parse(val, CultureInfo.InvariantCulture);
        });
        Console.WriteLine($"  preloaded {irtCache.Count} cached iRT predictions from {cachePath}");
    }
    TextWriter cacheLog = noIrt ? null : new StreamWriter(cachePath, append: true) { AutoFlush = false };

    Console.WriteLine($"--shardbuild: {allFastas.Length} databases -> {numShards} shards ({shardConcurrency}-concurrent " +
                      $"digest+write, {(noIrt ? "NO iRT" : "serialized iRT @ " + coreBudget + " threads, global cache")}). Out: {shardOutDir}");
    for (int i = 0; i < numShards; i++)
        Console.WriteLine($"  shard {i:D2}: {shardFastas[i].Count} dbs, {shardBytes[i] / 1048576} MB fasta");

    var shardSw = Stopwatch.StartNew();
    long[] shardEntryCounts = new long[numShards];
    bool[] shardSkipped = new bool[numShards];
    try
    {
        Parallel.For(0, numShards, new ParallelOptions { MaxDegreeOfParallelism = shardConcurrency }, si =>
        {
            string shardPath = Path.Combine(shardOutDir, $"bact_shard_{si:D2}.msl");
            // CHECKPOINT: a finalized .msl only exists on success (WriteStreaming merges spill sidecars at the
            // end), so its presence means this shard is done — skip it on a restart after an interruption.
            if (File.Exists(shardPath) && new FileInfo(shardPath).Length > 0)
            {
                shardSkipped[si] = true;
                Console.WriteLine($"  [shard {si:D2}] skip: already built ({new FileInfo(shardPath).Length / 1048576} MB)");
                return;
            }
            // Clear stale .spill~/.frags~ sidecars left by a prior crashed attempt at this shard.
            foreach (var tmp in Directory.GetFiles(shardOutDir, $"bact_shard_{si:D2}.msl.*~")) File.Delete(tmp);
            long cnt = 0;
            MslWriter.WriteStreaming(shardPath,
                ShardEntries(shardFastas[si], sharedPred, irtLock, coreBudget, irtCache, cacheLog, irtBusyTicks,
                    sFixedMods, sVariableMods, sDigest, () => cnt++));
            shardEntryCounts[si] = cnt;
            Console.WriteLine($"  [shard {si:D2}] done: {cnt} entries, {new FileInfo(shardPath).Length / 1048576} MB");
        });
    }
    finally { sharedPred?.Dispose(); cacheLog?.Dispose(); }
    shardSw.Stop();
    int skippedCount = shardSkipped.Count(b => b);
    double wallMin = shardSw.Elapsed.TotalMinutes;
    double irtBusyMin = TimeSpan.FromTicks(irtBusyTicks[0]).TotalMinutes;
    Console.WriteLine($"--shardbuild complete: {shardEntryCounts.Sum()} entries across {numShards - skippedCount} built shards" +
                      $"{(skippedCount > 0 ? $" ({skippedCount} skipped, already built)" : "")}, {irtCache.Count} unique iRT cached, " +
                      $"in {wallMin:F1} min.");
    if (!noIrt)
        Console.WriteLine($"  iRT predictor busy {irtBusyMin:F1} min of {wallMin:F1} min wall " +
                          $"({(wallMin > 0 ? 100 * (1 - irtBusyMin / wallMin) : 0):F0}% predictor-idle).");
    return 0;
}

if (positional.Length < 2)
{
    Console.Error.WriteLine("Usage: MslDbBuilder [--verify] <outputDir> <fasta1> [fasta2 ...]");
    return 1;
}

string outputDir = positional[0];
Directory.CreateDirectory(outputDir);
// FASTAs come from the remaining positionals, OR — to get past the OS command-line length limit when
// building thousands of databases — from a newline-separated list file passed as "@<path>".
string[] fastas;
if (positional.Length == 2 && positional[1].StartsWith("@"))
{
    fastas = File.ReadAllLines(positional[1].Substring(1))
        .Select(l => l.Trim()).Where(l => l.Length > 0).ToArray();
    Console.WriteLine($"Read {fastas.Length} FASTA paths from {positional[1].Substring(1)}");
}
else
{
    fastas = positional.Skip(1).ToArray();
}
Modification GetMod(string idWithMotif)
{
    var m = allMods.FirstOrDefault(x => x.IdWithMotif == idWithMotif);
    if (m == null) throw new Exception($"Modification '{idWithMotif}' not found in {ModsFile}");
    return m;
}
var fixedMods = new List<Modification> { GetMod("Carbamidomethyl on C") };
var variableMods = new List<Modification> { GetMod("Oxidation on M") };
var digestionParams = new DigestionParams(maxMissedCleavages: 0); // 0 MC: library defines the peptide space

Console.WriteLine($"Mods: fixed=[{string.Join(", ", fixedMods.Select(m => m.IdWithMotif))}] " +
                  $"variable=[{string.Join(", ", variableMods.Select(m => m.IdWithMotif))}]");
Console.WriteLine($"Converting {fastas.Length} database(s) -> {outputDir}\n");

// Chronologer predictor: store predicted iRT per entry so the .msl index is searchable by retention time.
// One instance reused across all databases (holds an unmanaged TorchSharp model — dispose at the end).
int rtThreads = Math.Max(1, Environment.ProcessorCount - 2);
using var rtPredictor = RetentionTimePredictorFactory.Create(PredictorType.Chronologer);
Console.WriteLine($"iRT predictor: {rtPredictor.GetType().Name} (maxThreads={rtThreads})");

int dbNum = 0;
foreach (var fasta in fastas)
{
    dbNum++;
    string name = Path.GetFileNameWithoutExtension(fasta);
    string outPath = Path.Combine(outputDir, name + ".msl");

    var proteins = ProteinDbLoader.LoadProteinFasta(
        fasta, generateTargets: true, decoyType: DecoyType.Reverse, isContaminant: false, out var errors);
    if (errors != null && errors.Count > 0)
        Console.Error.WriteLine($"  [{name}] {errors.Count} loader warning(s); first: {errors[0]}");

    int targetProteins = proteins.Count(p => !p.IsDecoy);
    int decoyProteins = proteins.Count(p => p.IsDecoy);

    // Dedup peptides by full sequence; TARGET WINS a target/decoy collision (drop the decoy).
    // Track the source-protein accession (for parsimony) and the peptide object (for iRT prediction).
    var byKey = new Dictionary<string, (LibrarySpectrum Spec, string Accession, PeptideWithSetModifications Pep)>();
    var products = new List<Product>();

    foreach (var protein in proteins)
    {
        foreach (var peptide in protein.Digest(digestionParams, fixedMods, variableMods))
        {
            // Skip peptides with an undefined mass (e.g. ambiguous residue 'X'): precursor m/z is NaN.
            if (!double.IsFinite(peptide.MonoisotopicMass)) continue;
            List<MatchedFragmentIon> ions;
            if (lean)
            {
                // Lean: store NO fragments. The search fragments this peptide on the fly.
                ions = new List<MatchedFragmentIon>();
            }
            else
            {
                products.Clear();
                peptide.Fragment(Dissociation, digestionParams.FragmentationTerminus, products);

                // Store fragment ions at charge 1 so the neutral mass is cleanly recoverable.
                ions = new List<MatchedFragmentIon>(products.Count);
                foreach (var p in products)
                {
                    if (double.IsNaN(p.NeutralMass)) continue;
                    ions.Add(new MatchedFragmentIon(p, p.NeutralMass.ToMz(1), 1.0, 1));
                }
                if (ions.Count == 0) continue;
            }

            double precursorMz = peptide.MonoisotopicMass.ToMz(NominalCharge);
            // rt is filled in below by a single batched iRT prediction over all kept peptides.
            var spec = new LibrarySpectrum(
                peptide.FullSequence, precursorMz, NominalCharge, ions, rt: null, isDecoy: protein.IsDecoy);

            if (!byKey.TryGetValue(peptide.FullSequence, out var existing))
                byKey[peptide.FullSequence] = (spec, protein.Accession, peptide);
            else if (existing.Spec.IsDecoy && !protein.IsDecoy)
                byKey[peptide.FullSequence] = (spec, protein.Accession, peptide); // target replaces a previously-seen decoy
            // otherwise keep existing (target stays; same-type: first wins)
        }
    }

    // Store predicted iRT (Chronologer) in each entry so the .msl index is searchable by retention time.
    // One batched inference over the whole de-duplicated set (targets AND decoys); assign back by sequence.
    var specBySeq = new Dictionary<string, LibrarySpectrum>(byKey.Count);
    foreach (var v in byKey.Values) specBySeq[v.Pep.FullSequence] = v.Spec;
    int irtFilled = 0;
    foreach (var r in rtPredictor.PredictRetentionTimeEquivalents(
                 byKey.Values.Select(v => (IRetentionPredictable)v.Pep), maxThreads: rtThreads))
    {
        if (r.PredictedValue is null) continue;
        if (specBySeq.TryGetValue(r.Peptide.FullSequence, out var sp)) { sp.RetentionTime = r.PredictedValue; irtFilled++; }
    }

    var spectra = byKey.Values.Select(v => v.Spec).ToList();
    int targetPep = spectra.Count(s => !s.IsDecoy);
    int decoyPep = spectra.Count(s => s.IsDecoy);

    // Convert to MslLibraryEntry and stamp the real source-protein accession, then write.
    var entries = new List<MslLibraryEntry>(byKey.Count);
    foreach (var (spec, accession, _) in byKey.Values)
    {
        var entry = MslLibraryEntry.FromLibrarySpectrum(spec);
        if (entry == null) continue;
        entry.ProteinAccession = accession;
        entries.Add(entry);
    }
    MslWriter.Write(outPath, entries);

    long bytes = new FileInfo(outPath).Length;
    Console.WriteLine($"[{dbNum}/{fastas.Length}] {name}.msl  " +
                      $"proteins {targetProteins}T/{decoyProteins}D  " +
                      $"peptides {targetPep}T/{decoyPep}D  entries {spectra.Count}  iRT {irtFilled}/{spectra.Count}  ({bytes / 1024} KB)");

    if (verify) VerifyRoundTrip(outPath, spectra);
}

Console.WriteLine($"\nDone. {fastas.Length} .msl libraries written to {outputDir}");
return 0;

// Reads the just-written .msl back through SpectralLibrary and compares every entry to what was
// written: deduped entry count, per-entry fragment count, fragment masses (max delta reveals
// float32-vs-double storage), and decoy flag.
static void VerifyRoundTrip(string path, List<LibrarySpectrum> written)
{
    using var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { path });
    var read = lib.GetAllLibrarySpectra().ToList();

    // MSL dedups by sequence/charge (last entry wins) — mirror that for an apples-to-apples count.
    var writtenByKey = new Dictionary<string, LibrarySpectrum>();
    foreach (var s in written) writtenByKey[s.Sequence + "/" + s.ChargeState] = s;

    int missing = 0, fragCountMismatch = 0, decoyMismatch = 0;
    double maxMassDelta = 0;
    foreach (var r in read)
    {
        if (!writtenByKey.TryGetValue(r.Sequence + "/" + r.ChargeState, out var w)) { missing++; continue; }
        if (r.IsDecoy != w.IsDecoy) decoyMismatch++;
        var rm = r.MatchedFragmentIons.Select(f => f.Mz).OrderBy(x => x).ToArray();
        var wm = w.MatchedFragmentIons.Select(f => f.Mz).OrderBy(x => x).ToArray();
        if (rm.Length != wm.Length) { fragCountMismatch++; continue; }
        for (int i = 0; i < rm.Length; i++)
        {
            double d = Math.Abs(rm[i] - wm[i]);
            if (d > maxMassDelta) maxMassDelta = d;
        }
    }

    bool ok = read.Count == writtenByKey.Count && missing == 0 && fragCountMismatch == 0 && decoyMismatch == 0;
    Console.WriteLine($"    verify: read {read.Count} vs written-deduped {writtenByKey.Count} | " +
                      $"missing {missing} | fragCountMismatch {fragCountMismatch} | decoyMismatch {decoyMismatch} | " +
                      $"maxFragMassDelta {maxMassDelta:E2} Da -> {(ok ? "ROUND-TRIP OK" : "MISMATCH")}");
}

// Reconstruct a searchable PeptideWithSetModifications from each .msl entry's FullSequence and
// re-fragment in DOUBLE precision. Confirms (a) every entry reconstructs, and (b) the recomputed
// double fragment COUNT matches the stored fragmentation — i.e. the library faithfully identifies
// the same peptides the search would fragment on the fly. This is the foundation of Strategy A.
// Parse confident target PSMs from a MetaMorpheus BasePSMs.psmtsv, predict iRT (SSRCalc3), regress
// observed RT vs predicted, and report the residual SD — the achievable RT-window width on this gradient.
static void RunRtCalib(string psmtsvPath, Dictionary<string, Modification> allKnownMods, PredictorType predType)
{
    using var reader = new StreamReader(psmtsvPath);
    string header = reader.ReadLine() ?? throw new Exception("empty psmtsv");
    var cols = header.Split('\t');
    int iFull = Array.IndexOf(cols, "Full Sequence");
    int iRt = Array.IndexOf(cols, "Scan Retention Time");
    int iQ = Array.IndexOf(cols, "QValue");
    int iDct = Array.IndexOf(cols, "Decoy/Contaminant/Target");
    if (iFull < 0 || iRt < 0 || iQ < 0 || iDct < 0) throw new Exception("missing required columns");

    var predictor = RetentionTimePredictorFactory.Create(predType);
    var obs = new List<double>();
    var pred = new List<double>();
    int rows = 0, predFailed = 0, parseFailed = 0;

    string? line;
    while ((line = reader.ReadLine()) != null)
    {
        var c = line.Split('\t');
        if (c.Length <= Math.Max(iFull, Math.Max(iRt, Math.Max(iQ, iDct)))) continue;
        if (c[iDct] != "T") continue;
        if (!double.TryParse(c[iQ], out double q) || q >= 0.01) continue;
        if (!double.TryParse(c[iRt], out double rt)) continue;
        string full = c[iFull];
        if (full.Contains('|')) continue; // ambiguous
        rows++;

        Proteomics.ProteolyticDigestion.PeptideWithSetModifications pep;
        try { pep = new Proteomics.ProteolyticDigestion.PeptideWithSetModifications(full, allKnownMods); }
        catch { parseFailed++; continue; }

        double? p = predictor.PredictRetentionTimeEquivalent(pep, out _);
        if (p is null) { predFailed++; continue; }
        obs.Add(rt); pred.Add(p.Value);
    }

    int n = obs.Count;
    if (n < 2) { Console.WriteLine($"RTCALIB: too few usable PSMs (n={n})"); return; }

    // OLS: observed = slope*predicted + intercept
    double mx = pred.Average(), my = obs.Average();
    double sxy = 0, sxx = 0, syy = 0;
    for (int i = 0; i < n; i++)
    {
        double dx = pred[i] - mx, dy = obs[i] - my;
        sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    double slope = sxy / sxx, intercept = my - slope * mx;
    double r2 = (sxy * sxy) / (sxx * syy);

    var resid = new List<double>(n);
    for (int i = 0; i < n; i++) resid.Add(obs[i] - (slope * pred[i] + intercept));
    double residSd = Math.Sqrt(resid.Sum(r => r * r) / n);
    var absResSorted = resid.Select(Math.Abs).OrderBy(x => x).ToList();
    double P(double frac) => absResSorted[(int)(frac * (n - 1))];
    double rtMin = obs.Min(), rtMax = obs.Max(), span = rtMax - rtMin;

    Console.WriteLine($"RTCALIB ({predType}): n={n} parseFail={parseFailed} predFail={predFailed}");
    Console.WriteLine($"  OLS observedRT = {slope:F4}*iRT + {intercept:F2}   R2={r2:F4}");
    Console.WriteLine($"  residual SD = {residSd:F2} min   |resid| p50={P(0.5):F2} p90={P(0.9):F2} p95={P(0.95):F2} p99={P(0.99):F2}");
    Console.WriteLine($"  gradient span = {span:F1} min ({rtMin:F1}-{rtMax:F1})");
    double win95 = 2 * P(0.95);
    Console.WriteLine($"  => RT window for 95% containment = +/-{P(0.95):F2} min = {win95:F2}/{span:F1} = {100.0 * win95 / span:F1}% of gradient (selectivity x{span / win95:F1})");
}

static void ReadAndCheckPeptides(string mslPath, Dictionary<string, Modification> allKnownMods)
{
    using var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { mslPath });
    var products = new List<Product>();

    int entries = 0, reconstructed = 0, parseFailed = 0, fragCountMismatch = 0, decoys = 0;
    long totalDoubleFragments = 0;
    int samplesShown = 0;

    foreach (var spec in lib.GetAllLibrarySpectra())
    {
        entries++;
        if (spec.IsDecoy) decoys++;

        Proteomics.ProteolyticDigestion.PeptideWithSetModifications pep;
        try { pep = new Proteomics.ProteolyticDigestion.PeptideWithSetModifications(spec.Sequence, allKnownMods); }
        catch { parseFailed++; continue; }
        reconstructed++;

        products.Clear();
        pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
        int dbl = products.Count(p => !double.IsNaN(p.NeutralMass));
        totalDoubleFragments += dbl;

        // The stored ions were generated by the same HCD/both-termini fragmentation, so the
        // recomputed double-fragment count should equal the stored ion count.
        if (dbl != spec.MatchedFragmentIons.Count) fragCountMismatch++;

        if (samplesShown < 3)
        {
            Console.WriteLine($"    e.g. {spec.Sequence} | stored {spec.MatchedFragmentIons.Count} ions, " +
                              $"recomputed {dbl} double fragments | decoy={spec.IsDecoy}");
            samplesShown++;
        }
    }

    bool ok = parseFailed == 0 && fragCountMismatch == 0 && reconstructed == entries;
    Console.WriteLine($"{Path.GetFileName(mslPath)}: entries {entries} ({decoys} decoy) | reconstructed {reconstructed} | " +
                      $"parseFailed {parseFailed} | fragCountMismatch {fragCountMismatch} | " +
                      $"totalDoubleFragments {totalDoubleFragments} -> {(ok ? "READ+REFRAGMENT OK" : "ISSUES")}");
}

// Lazily yields one shard's LEAN entries, one database at a time, so WriteStreaming never holds them all
// in RAM. Per db: load (target+decoy), digest, dedup by sequence (target wins), batched iRT, yield each
// peptide as a fragment-less MslLibraryEntry tagged "db|accession".
static IEnumerable<MslLibraryEntry> ShardEntries(
    List<string> fastas, IRetentionTimePredictor predictor, object irtLock, int irtThreads,
    ShardedSeqCache irtCache, TextWriter cacheLog, long[] irtBusyTicks,
    List<Modification> fixedMods, List<Modification> variableMods, DigestionParams dp, Action onEntry)
{
    foreach (var fasta in fastas)
    {
        string db = Path.GetFileNameWithoutExtension(fasta);
        var proteins = ProteinDbLoader.LoadProteinFasta(fasta, generateTargets: true,
            decoyType: DecoyType.Reverse, isContaminant: false, out _);

        // Dedup by sequence within this proteome, but KEEP EVERY source accession (a peptide shared by
        // several proteins of the same organism must list them all, else transient protein parsimony is
        // wrong). Cross-proteome duplicates are NOT collapsed — each db emits its own entry (the search
        // unions per "db|acc" group), so organism attribution is preserved. Target wins over decoy on a
        // sequence collision (and resets the accession list to the target's proteins).
        var byKey = new Dictionary<string, (List<string> Accs, PeptideWithSetModifications Pep, bool Decoy)>();
        foreach (var protein in proteins)
            foreach (var pep in protein.Digest(dp, fixedMods, variableMods))
            {
                if (!byKey.TryGetValue(pep.FullSequence, out var ex))
                    byKey[pep.FullSequence] = (new List<string> { protein.Accession }, pep, protein.IsDecoy);
                else if (ex.Decoy && !protein.IsDecoy)
                    byKey[pep.FullSequence] = (new List<string> { protein.Accession }, pep, false); // target supersedes decoys
                else if (ex.Decoy == protein.IsDecoy && !ex.Accs.Contains(protein.Accession))
                    ex.Accs.Add(protein.Accession); // same class: accumulate the source protein (mutates the shared list)
            }

        // LAZY GLOBAL iRT: predict only sequences not seen in ANY prior database, batched under the shared
        // lock (the one predictor saturates all cores; predictor == null => --noirt build, rt stays null).
        if (predictor != null)
        {
            var missSeqs = new List<string>();
            foreach (var seq in byKey.Keys)
                if (!irtCache.ContainsKey(seq)) missSeqs.Add(seq);
            if (missSeqs.Count > 0)
            {
                lock (irtLock)
                {
                    // Re-check under the lock: another shard may have predicted some of these while we waited.
                    var toPredict = new List<PeptideWithSetModifications>(missSeqs.Count);
                    foreach (var seq in missSeqs)
                        if (!irtCache.ContainsKey(seq)) toPredict.Add(byKey[seq].Pep);
                    if (toPredict.Count > 0)
                    {
                        var fresh = new List<KeyValuePair<string, double?>>(toPredict.Count);
                        var returned = new HashSet<string>();
                        var psw = Stopwatch.StartNew();
                        foreach (var r in predictor.PredictRetentionTimeEquivalents(
                                     toPredict.Select(p => (IRetentionPredictable)p), maxThreads: irtThreads))
                        {
                            irtCache[r.Peptide.FullSequence] = r.PredictedValue;
                            returned.Add(r.Peptide.FullSequence);
                            fresh.Add(new KeyValuePair<string, double?>(r.Peptide.FullSequence, r.PredictedValue));
                        }
                        psw.Stop();
                        Interlocked.Add(ref irtBusyTicks[0], psw.ElapsedTicks);
                        // Cache sequences the predictor returned nothing for as null, so they aren't re-attempted.
                        foreach (var p in toPredict)
                            if (returned.Add(p.FullSequence))
                            {
                                irtCache[p.FullSequence] = null;
                                fresh.Add(new KeyValuePair<string, double?>(p.FullSequence, null));
                            }
                        if (cacheLog != null)
                        {
                            foreach (var kv in fresh)
                                cacheLog.WriteLine(kv.Key + "\t" + (kv.Value?.ToString("R", CultureInfo.InvariantCulture) ?? ""));
                            cacheLog.Flush();
                        }
                    }
                }
            }
        }

        foreach (var kv in byKey)
        {
            var (accs, pep, decoy) = kv.Value;
            // Skip peptides with an undefined mass (e.g. sequences containing the ambiguous residue 'X'):
            // their precursor m/z is NaN, which is useless for search and pollutes the m/z index.
            if (!double.IsFinite(pep.MonoisotopicMass)) continue;
            double? rt = irtCache.TryGetValue(kv.Key, out var v) ? v : null;
            // LEAN entries carry NO fragments (the search fragments candidates on the fly), so build the
            // MslLibraryEntry POCO directly and skip LibrarySpectrum/MzSpectrum — that ctor does two LINQ
            // Select().ToArray() over the peaks, an Array.Sort, and allocates an empty peaks List, all pure
            // waste here and multiplied by billions of entries. Field mapping mirrors FromLibrarySpectrum
            // (BaseSequence via pep, already computed); the rest keep their POCO defaults (empty fragment
            // list, IonMobility 0, QValue NaN, ElutionGroupId 0, Peptide/Predicted/Unknown/Nce 0).
            var entry = new MslLibraryEntry
            {
                FullSequence = pep.FullSequence,
                BaseSequence = pep.BaseSequence,
                PrecursorMz = pep.MonoisotopicMass.ToMz(2),
                ChargeState = 2,
                RetentionTime = rt ?? 0.0,
                IsDecoy = decoy,
                ProteinAccession = db + "|" + string.Join(";", accs),
                ProteinName = string.Empty,
                GeneName = string.Empty,
            };
            onEntry();
            yield return entry;
        }
    }
}

// A single ConcurrentDictionary cannot reliably hold the full bacterial unique-sequence set (~1-2 billion
// entries): its internal bucket array hits a .NET array/scale limit and corrupts on resize (NRE in
// TryAddInternal) even with hundreds of GB of RAM free. Shard by key hash into N independent dictionaries so
// each holds only ~total/N entries (tens of millions — well within limits) and every resize is a small
// allocation. N is a power of two so the shard index is a cheap mask. Same lock-free concurrent-read /
// concurrent-write semantics as the underlying ConcurrentDictionary, per shard.
sealed class ShardedSeqCache
{
    private readonly System.Collections.Concurrent.ConcurrentDictionary<string, double?>[] _shards;
    private readonly int _mask;

    public ShardedSeqCache(int shardCount = 128)
    {
        _shards = new System.Collections.Concurrent.ConcurrentDictionary<string, double?>[shardCount];
        for (int i = 0; i < shardCount; i++)
            _shards[i] = new System.Collections.Concurrent.ConcurrentDictionary<string, double?>();
        _mask = shardCount - 1;   // shardCount must be a power of two
    }

    private System.Collections.Concurrent.ConcurrentDictionary<string, double?> ShardFor(string key)
        => _shards[(key.GetHashCode() & 0x7fffffff) & _mask];

    public bool ContainsKey(string key) => ShardFor(key).ContainsKey(key);
    public bool TryGetValue(string key, out double? value) => ShardFor(key).TryGetValue(key, out value);
    public double? this[string key] { set => ShardFor(key)[key] = value; }
    public long Count { get { long c = 0; foreach (var s in _shards) c += s.Count; return c; } }
}
