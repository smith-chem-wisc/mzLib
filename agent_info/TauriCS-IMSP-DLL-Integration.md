# TauriCS IMSP DLL Integration Guide

This guide covers using the mzLib IMSP export services from a TauriCS app where the backend is C# and the frontend is TypeScript/Rust/Tauri.

The recommended integration is:

1. Build or publish mzLib DLLs.
2. Reference those DLLs from the C# backend project.
3. Expose small C# backend commands for `mzML -> .imsp file` and `mzML -> IMSP bytes`.
4. Call those backend commands from the TypeScript frontend.

Do not try to load `MassSpectrometry.dll` or `Readers.dll` directly from Rust with `dlopen`/FFI. These are managed .NET assemblies, not C ABI native libraries. Rust should talk to the C# backend layer, or to a separate .NET sidecar process, unless you build a dedicated NativeAOT/C ABI wrapper.

## DLLs Needed

For the path-based mzML workflow, reference these assemblies:

- `Readers.dll`
- `MassSpectrometry.dll`
- `MzLibUtil.dll`
- `Chemistry.dll`
- `Transcriptomics.dll`
- `UsefulProteomicsDatabases.dll`

`Readers.dll` contains the mzML file-path adapter:

```csharp
Readers.MzmlImspExportService
```

`MassSpectrometry.dll` contains the core IMSP writer:

```csharp
MassSpectrometry.ImspExportService
```

If your backend already has `MsDataScan` objects and does not need to load mzML files, you can use only the core `MassSpectrometry` layer and its dependencies.

## Build DLLs

From the repository root:

```powershell
dotnet publish .\mzLib\Readers\Readers.csproj -c Release -r win-x64 --self-contained false
```

For Intel Mac:

```powershell
dotnet publish .\mzLib\Readers\Readers.csproj -c Release -r osx-x64 --self-contained false
```

The DLLs will be in:

```text
mzLib/Readers/bin/Release/net8.0/<runtime>/publish/
```

For example:

```text
mzLib/Readers/bin/Release/net8.0/win-x64/publish/
mzLib/Readers/bin/Release/net8.0/osx-x64/publish/
```

Use the publish directory as the source for DLL references or copy the needed DLLs into a stable `lib/mzLib` folder in your TauriCS backend project.

## Reference DLLs From The C# Backend

In the C# backend `.csproj`, reference the DLLs with `HintPath`.

Example:

```xml
<ItemGroup>
  <Reference Include="Readers">
    <HintPath>lib/mzLib/Readers.dll</HintPath>
    <Private>true</Private>
  </Reference>
  <Reference Include="MassSpectrometry">
    <HintPath>lib/mzLib/MassSpectrometry.dll</HintPath>
    <Private>true</Private>
  </Reference>
  <Reference Include="MzLibUtil">
    <HintPath>lib/mzLib/MzLibUtil.dll</HintPath>
    <Private>true</Private>
  </Reference>
  <Reference Include="Chemistry">
    <HintPath>lib/mzLib/Chemistry.dll</HintPath>
    <Private>true</Private>
  </Reference>
  <Reference Include="Transcriptomics">
    <HintPath>lib/mzLib/Transcriptomics.dll</HintPath>
    <Private>true</Private>
  </Reference>
  <Reference Include="UsefulProteomicsDatabases">
    <HintPath>lib/mzLib/UsefulProteomicsDatabases.dll</HintPath>
    <Private>true</Private>
  </Reference>
</ItemGroup>
```

If the TauriCS backend can reference projects directly during development, a project reference is cleaner:

```xml
<ItemGroup>
  <ProjectReference Include="..\..\mzLib\Readers\Readers.csproj" />
</ItemGroup>
```

For app packaging, published DLLs are usually easier to control.

## C# Backend Service

Create a small app-level wrapper around mzLib. This keeps mzLib-specific details out of your Tauri command layer.

```csharp
using Readers;

public sealed class ImspBackendService
{
    private readonly IMzmlImspExportService imspExportService;

    public ImspBackendService()
        : this(new MzmlImspExportService())
    {
    }

    public ImspBackendService(IMzmlImspExportService imspExportService)
    {
        this.imspExportService = imspExportService;
    }

    public string ConvertMzmlToImspFile(string mzmlPath, string? outputPath = null,
        int binsPerDalton = MassSpectrometry.ImspExportService.DefaultBinsPerDalton,
        double intensityThreshold = MassSpectrometry.ImspExportService.DefaultIntensityThreshold)
    {
        return imspExportService.ConvertToImspFile(
            mzmlPath,
            outputPath,
            binsPerDalton,
            intensityThreshold);
    }

    public byte[] ConvertMzmlToImspBytes(string mzmlPath,
        int binsPerDalton = MassSpectrometry.ImspExportService.DefaultBinsPerDalton,
        double intensityThreshold = MassSpectrometry.ImspExportService.DefaultIntensityThreshold)
    {
        return imspExportService.ConvertToImspBytes(
            mzmlPath,
            binsPerDalton,
            intensityThreshold);
    }
}
```

## Expose TauriCS Commands

The exact attribute and command shape depends on the TauriCS setup in your app. The command layer should be thin and should call the backend service above.

Example shape:

```csharp
public sealed class ImspCommands
{
    private readonly ImspBackendService service = new();

    public string ConvertMzmlToImspFile(string mzmlPath, string? outputPath = null,
        int binsPerDalton = 100, double intensityThreshold = 10000)
    {
        return service.ConvertMzmlToImspFile(
            mzmlPath,
            outputPath,
            binsPerDalton,
            intensityThreshold);
    }

    public string ConvertMzmlToImspBase64(string mzmlPath,
        int binsPerDalton = 100, double intensityThreshold = 10000)
    {
        byte[] bytes = service.ConvertMzmlToImspBytes(
            mzmlPath,
            binsPerDalton,
            intensityThreshold);

        return Convert.ToBase64String(bytes);
    }
}
```

Returning a file path is preferred for large files. Returning base64 is convenient for small tests, but it increases payload size and memory use. For production-sized IMSP files, write the `.imsp` to disk and return the path.

## TypeScript Frontend Calls

The TypeScript call depends on how TauriCS maps C# commands into Tauri. The shape is typically similar to Tauri `invoke`.

File output example:

```ts
const imspPath = await invoke<string>("convert_mzml_to_imsp_file", {
  mzmlPath: "C:\\data\\sample.mzML",
  outputPath: "C:\\data\\sample.imsp",
  binsPerDalton: 100,
  intensityThreshold: 10000,
});
```

Byte output example:

```ts
const imspBase64 = await invoke<string>("convert_mzml_to_imsp_base64", {
  mzmlPath: "C:\\data\\sample.mzML",
  binsPerDalton: 100,
  intensityThreshold: 10000,
});

const bytes = Uint8Array.from(atob(imspBase64), c => c.charCodeAt(0));
```

For large files, avoid base64 and use the file path workflow.

## Rust/Tauri Notes

Use Rust for shell integration, app lifecycle, file dialogs, and command routing. Let the C# backend do mzML reading and IMSP writing.

If a Rust command needs to trigger IMSP export, route the request to the C# backend rather than loading mzLib DLLs directly.

If you eventually need a pure Rust-to-library boundary, there are two realistic options:

- Run a small .NET sidecar executable and communicate with stdin/stdout, HTTP, named pipes, or local sockets.
- Build a dedicated NativeAOT wrapper that exports a C ABI. That is a separate project and requires careful packaging/testing.

## Mac Compatibility

For Intel Mac, publish the DLL set with:

```powershell
dotnet publish .\mzLib\Readers\Readers.csproj -c Release -r osx-x64 --self-contained false
```

The mzML-only path is the safest cross-platform target. Thermo RAW and vendor-specific readers depend on native/vendor libraries and should not be assumed to work on macOS.

If packaging for macOS, test with a real `.mzML` file on the target machine:

```csharp
var service = new Readers.MzmlImspExportService();
var imspPath = service.ConvertToImspFile("/Users/alex/data/sample.mzML");
```

## Error Handling

The current mzML adapter throws exceptions for invalid inputs:

- `ArgumentException` for empty paths.
- `NotSupportedException` for non-`.mzML` inputs.
- File/reader exceptions if the mzML cannot be opened or parsed.

Wrap command bodies in `try/catch` and return structured errors to the frontend.

Example:

```csharp
try
{
    return service.ConvertMzmlToImspFile(mzmlPath, outputPath);
}
catch (Exception ex)
{
    throw new InvalidOperationException($"IMSP export failed: {ex.Message}", ex);
}
```

## Recommended First Integration Test

1. Add the DLL references to the C# backend.
2. Call `new MzmlImspExportService().ConvertToImspFile(...)` from a backend unit test or temporary command.
3. Verify the returned path exists and starts with the IMSP magic bytes.

```csharp
byte[] header = File.ReadAllBytes(imspPath).Take(4).ToArray();
bool isImsp = header.SequenceEqual(new byte[] { (byte)'I', (byte)'M', (byte)'S', (byte)'P' });
```

Once that works, wire the same backend call to the TypeScript UI.
