#r "C:\\Users\\Nic\\source\\repos\\mzLib\\mzLib\\Omics\\bin\\Debug\\net8.0\\Omics.dll"
#r "C:\\Users\\Nic\\source\\repos\\mzLib\\mzLib\\Chemistry\\bin\\Debug\\net8.0\\Chemistry.dll"
using System;
using Omics.Modifications;

var mod = Mods.GetModification("Water Loss on D", ModificationNamingConvention.Mixed);
Console.WriteLine(mod?.IdWithMotif ?? "null IdWithMotif");
Console.WriteLine(mod?.OriginalId ?? "null OriginalId");
