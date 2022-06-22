using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test
{
    class aFindingSpecificMassesInScans
    {
        [Test]
        // Never finished as I just did this in the other program I was working on. Hopefully this doesnt come back to bite me D:
        public void FindSpecificMasses() 
        {
            // Loads in proteins from a scan and creates a List of doubles representing their m/z values with charges up to 60 from their monoisotopic mass
            string filename = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\Six_Protein_Standard.fasta";
            var proteinList = ProteinDbLoader.LoadProteinFasta(filename, true, DecoyType.None, false, out var dbErrors);
            List<Peptide> proteinsWithMasses = new();
            foreach (var protein in proteinList)
            {
                Peptide prot = new(protein.BaseSequence);
                proteinsWithMasses.Add(prot);
            }
            int[] charges = new int[60];
            for (int i = 1; i < 60; i++)
            {
                charges[i-1] = i;
            }
            double[] masses = proteinsWithMasses.Select(p => p.MonoisotopicMass).ToArray();
            List<double> ions = new();
            for (int j = 0; j < masses.Count(); j++)
            {
                for (int i = 0; i < charges.Length; i++)
                {
                    ions.Add(masses[j] / charges[i]);
                }
                ions.Remove(double.PositiveInfinity);
            }

            int breakpoint = 0;

            //Load in masses
            filename = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\2021-01-06_TopDownStandard_YeastFraction1.raw";
            List<MsDataScan> scans = ThermoRawFileReader.LoadAllStaticData(filename).GetAllScansList();
            foreach (var scan in scans)
            {

            }

        }
    }
}
