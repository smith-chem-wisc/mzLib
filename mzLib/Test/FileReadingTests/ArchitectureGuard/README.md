# Architecture guard fixture

`ThermoFisher.CommonCore.MassPrecisionEstimator.dll` here is **test data, not a dependency**. Nothing
references it, it is not shipped in the NuGet package, and it is not the copy the Thermo readers use.

It exists because `TestNoArchitectureSpecificAssemblyReferences` scans the build output for references to
assemblies stamped for a single architecture, and a scan like that passes just as happily when its logic
is broken as when there is nothing to find. `TheClassifierRecognisesAnArchitectureStampedAssembly` uses
this file as the positive control: it is IL-only with its PE machine field stamped `AMD64`, so the
classifier must report it as architecture-specific.

**This copy is pinned deliberately, rather than read from `Readers/Thermo/`.** If the vendored Thermo
binaries are ever updated — including to an architecture-neutral build, which would be good news — the
control still has a stamped assembly to assert against, and nobody has to choose between updating the
package and keeping the test honest.

So: do not "fix" this file. If it ever stops being `AMD64`-stamped, the control is measuring nothing, and
replacing it with any other architecture-stamped managed assembly is the repair.
