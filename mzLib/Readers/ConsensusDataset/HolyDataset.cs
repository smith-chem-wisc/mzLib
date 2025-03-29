using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Media;
using System.Text;
using System.Threading.Tasks;
using Readers;
using System.IO;
using TopDownProteomics.IO.Resid;

namespace Readers.ConsensusDataset
{
    public class HolyDatasetMST
    {
        /*
         * Returns output given by outPath, user Defined
         */
        public HolyDatasetMST(string exePath, string spectraPath, string dataPath,
            string outPath) //TODO given this has never ran to completion, more tests needed.
        {
            exePath = @"" + exePath;
            spectraPath = @"" + spectraPath;
            dataPath = @"" + dataPath;
            outPath = @"" + outPath;
            string exeParams = " -s " + "\"" + spectraPath + "\"" + " -d " + "\"" + dataPath + "\"" + " -o " + "\"" +
                               outPath + "\"" +
                               " -ic 2 -f 20 -MinLength 7 -MaxLength 1000000 -MinCharge 1 -MaxCharge 60 -MinFragCharge 1 -MaxFragCharge 10 -MinMass 0 -MaxMass 30000 -tda 1"; //100000
            var process = new Process
            {
                StartInfo =
                {
                    FileName = exePath,
                    WorkingDirectory = Path.GetDirectoryName(exePath),
                    Arguments = exeParams,
                }

            };

            process.Start();

            process.WaitForExit();


        }

        /*
         * OutPath left empty, use for result handling purposes
         */
        public HolyDatasetMST(string exePath, string spectraPath, string dataPath) : this(exePath, spectraPath,
            dataPath,
            Path.GetDirectoryName(@"" + spectraPath))
        {
            //result handling
            MsPathFinderTResultFile
                result = new MsPathFinderTResultFile(Path.ChangeExtension(spectraPath,
                    "_IcTda.tsv")); //suppose spectraPath.raw -> spectraPath_IcTda.tsv. TODO double check. Small thing.
            var dict = result.ToDictListList();


        }



        class userSimulation //what would the user have in main? TODO remove, for reference only.
        {
            string exePath = @"C:\Program Files\Informed-Proteomics\MSPathFinderT";
            string spectraPath = @"C:\Users\avnib\Desktop\SEOutput\RAW\Ecoli_SEC4_F6.raw";

            string datPath =
                @"C:\Users\avnib\Desktop\Databases\uniprotkb_proteome_UP000005640_2025_03_11.fasta\uniprotkb_proteome_UP000005640_2025_03_11.fasta";

            string output = @"C:\Users\avnib\Desktop\SEOutput\MST";
            string tsv = @"C:\Users\avnib\Desktop\02-17-20_jurkat_td_rep1_fract1_IcTda.tsv";


        }



    }
}
