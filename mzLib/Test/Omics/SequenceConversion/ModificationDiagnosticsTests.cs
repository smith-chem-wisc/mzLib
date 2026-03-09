using System;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;

namespace Test.Omics.SequenceConversion;

/// <summary>
/// Diagnostic tests to understand modification data structure
/// </summary>
[TestFixture]
public class ModificationDiagnosticsTests
{
    [Test]
    public void DiagnoseMods_PrintSampleModifications()
    {
        // Print first 10 protein mods
        Console.WriteLine("=== SAMPLE PROTEIN MODIFICATIONS ===");
        foreach (var mod in Mods.AllProteinModsList.Take(10))
        {
            Console.WriteLine($"IdWithMotif: {mod.IdWithMotif}");
            Console.WriteLine($"OriginalId: {mod.OriginalId}");
            Console.WriteLine($"ModificationType: {mod.ModificationType}");
            Console.WriteLine($"Accession: {mod.Accession}");
            if (mod.DatabaseReference != null && mod.DatabaseReference.Count > 0)
            {
                foreach (var kvp in mod.DatabaseReference)
                {
                    Console.WriteLine($"  DatabaseReference[{kvp.Key}]: {string.Join(", ", kvp.Value)}");
                }
            }
            Console.WriteLine();
        }

        // Print first 10 RNA mods
        Console.WriteLine("\n=== SAMPLE RNA MODIFICATIONS ===");
        foreach (var mod in Mods.AllRnaModsList.Take(10))
        {
            Console.WriteLine($"IdWithMotif: {mod.IdWithMotif}");
            Console.WriteLine($"OriginalId: {mod.OriginalId}");
            Console.WriteLine($"ModificationType: {mod.ModificationType}");
            Console.WriteLine($"Accession: {mod.Accession}");
            if (mod.DatabaseReference != null && mod.DatabaseReference.Count > 0)
            {
                foreach (var kvp in mod.DatabaseReference)
                {
                    Console.WriteLine($"  DatabaseReference[{kvp.Key}]: {string.Join(", ", kvp.Value)}");
                }
            }
            Console.WriteLine();
        }
    }

    [Test]
    public void DiagnoseMods_FindOxidationMod()
    {
        Console.WriteLine("=== SEARCHING FOR OXIDATION MODS ===");
        
        var oxidationMods = Mods.AllKnownMods
            .Where(m => m.IdWithMotif?.Contains("Oxidation", StringComparison.OrdinalIgnoreCase) == true
                     || m.OriginalId?.Contains("Oxidation", StringComparison.OrdinalIgnoreCase) == true)
            .ToList();

        Console.WriteLine($"Found {oxidationMods.Count} oxidation mods");
        
        foreach (var mod in oxidationMods.Take(5))
        {
            Console.WriteLine($"\nIdWithMotif: {mod.IdWithMotif}");
            Console.WriteLine($"OriginalId: {mod.OriginalId}");
            Console.WriteLine($"ModificationType: {mod.ModificationType}");
            Console.WriteLine($"Accession: {mod.Accession}");
            if (mod.DatabaseReference != null && mod.DatabaseReference.Count > 0)
            {
                foreach (var kvp in mod.DatabaseReference)
                {
                    Console.WriteLine($"  DatabaseReference[{kvp.Key}]: {string.Join(", ", kvp.Value)}");
                }
            }
        }
    }

    [Test]
    public void DiagnoseMods_FindUnimodReference()
    {
        Console.WriteLine("=== SEARCHING FOR MODS WITH UNIMOD REFERENCES ===");
        
        var unimodMods = Mods.AllKnownMods
            .Where(m => m.DatabaseReference?.ContainsKey("UNIMOD") == true)
            .Take(10)
            .ToList();

        Console.WriteLine($"Found {unimodMods.Count} mods with UNIMOD references (showing first 10)");
        
        foreach (var mod in unimodMods)
        {
            Console.WriteLine($"\nIdWithMotif: {mod.IdWithMotif}");
            Console.WriteLine($"OriginalId: {mod.OriginalId}");
            Console.WriteLine($"ModificationType: {mod.ModificationType}");
            if (mod.DatabaseReference.TryGetValue("UNIMOD", out var unimodIds))
            {
                Console.WriteLine($"  UNIMOD IDs: {string.Join(", ", unimodIds)}");
            }
        }
    }

    [Test]
    public void DiagnoseMods_FindTerminalMods()
    {
        Console.WriteLine("=== SEARCHING FOR TERMINAL MODIFICATIONS ===");
        
        var terminalMods = Mods.AllKnownMods
            .Where(m => m.LocationRestriction != null && 
                       (m.LocationRestriction.Contains("N-terminal") || 
                        m.LocationRestriction.Contains("C-terminal")))
            .Take(10)
            .ToList();

        Console.WriteLine($"Found {terminalMods.Count} terminal mods (showing first 10)");
        
        foreach (var mod in terminalMods)
        {
            Console.WriteLine($"\nIdWithMotif: {mod.IdWithMotif}");
            Console.WriteLine($"OriginalId: {mod.OriginalId}");
            Console.WriteLine($"LocationRestriction: {mod.LocationRestriction}");
            Console.WriteLine($"Target: {mod.Target}");
        }
    }

    [Test]
    public void DiagnoseMods_TestGetModificationMethods()
    {
        Console.WriteLine("=== TESTING GetModification METHODS ===");

        // Test by IdWithMotif
        var oxidation = Mods.GetModification("Oxidation on M");
        Console.WriteLine($"GetModification('Oxidation on M'): {oxidation?.IdWithMotif ?? "NULL"}");

        // Test by full ID with type prefix
        var oxidation2 = Mods.GetModification("Common Variable:Oxidation on M");
        Console.WriteLine($"GetModification('Common Variable:Oxidation on M'): {oxidation2?.IdWithMotif ?? "NULL"}");

        // Test with Mixed convention
        var oxidation3 = Mods.GetModification("Oxidation on M", ModificationNamingConvention.Mixed);
        Console.WriteLine($"GetModification('Oxidation on M', Mixed): {oxidation3?.IdWithMotif ?? "NULL"}");

        // Test UNIMOD ID
        var unimod21 = Mods.GetModification("UNIMOD:21", ModificationNamingConvention.Mixed);
        Console.WriteLine($"GetModification('UNIMOD:21', Mixed): {unimod21?.IdWithMotif ?? "NULL"}");

        // Test m6A RNA mod
        Console.WriteLine($"\n=== SEARCHING FOR m6A RNA MOD ===");
        var m6A = Mods.GetModification("m6A", ModificationNamingConvention.Mixed);
        Console.WriteLine($"GetModification('m6A', Mixed): {m6A?.IdWithMotif ?? "NULL"}");
        
        // Search for m6A in all mods
        var m6AMods = Mods.AllKnownMods
            .Where(m => m.IdWithMotif?.Contains("m6A") == true || m.OriginalId?.Contains("m6A") == true)
            .ToList();
        Console.WriteLine($"Mods containing 'm6A': {m6AMods.Count}");
        foreach (var mod in m6AMods.Take(5))
        {
            Console.WriteLine($"  - IdWithMotif: {mod.IdWithMotif}, OriginalId: {mod.OriginalId}, ModificationType: {mod.ModificationType}");
        }

        // Try to find phosphorylation (UNIMOD:21)
        var phospho = Mods.AllKnownMods.FirstOrDefault(m => 
            m.DatabaseReference?.ContainsKey("Unimod") == true &&
            m.DatabaseReference["Unimod"].Contains("21"));
        Console.WriteLine($"\nMod with UNIMOD:21 reference:");
        if (phospho != null)
        {
            Console.WriteLine($"  IdWithMotif: {phospho.IdWithMotif}");
            Console.WriteLine($"  OriginalId: {phospho.OriginalId}");
            Console.WriteLine($"  ModificationType: {phospho.ModificationType}");
        }
        else
        {
            Console.WriteLine("  NOT FOUND");
        }
    }
}
