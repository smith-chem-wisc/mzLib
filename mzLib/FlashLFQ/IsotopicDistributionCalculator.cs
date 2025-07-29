using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public sealed class IsotopicDistributionCalculator
    {
        internal static double AveragineMass; // Averagine, for protEINS
        internal static double AveratideMass; // Averatide, for ribonucleoTIDES (RNA)

        // Magic numbers determined by https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/
        // These number are for peptides/proteins
        internal static double AverageC = 4.9384;
        internal static double AverageH = 7.7583;
        internal static double AverageO = 1.4773;
        internal static double AverageN = 1.3577;
        internal static double AverageS = 0.0417;

        // These numbers are for RNA
        internal static double RnaAverageC = 9.5;
        internal static double RnaAverageH = 12.75;
        internal static double RnaAverageO = 3.75;
        internal static double RnaAverageN = 5;

        static IsotopicDistributionCalculator()
        {
            AveragineMass =
                PeriodicTable.GetElement("C").AverageMass * AverageC +
                PeriodicTable.GetElement("H").AverageMass * AverageH +
                PeriodicTable.GetElement("O").AverageMass * AverageO +
                PeriodicTable.GetElement("N").AverageMass * AverageN +
                PeriodicTable.GetElement("S").AverageMass * AverageS;

            AveratideMass =
                PeriodicTable.GetElement("C").AverageMass * RnaAverageC +
                PeriodicTable.GetElement("H").AverageMass * RnaAverageH +
                PeriodicTable.GetElement("O").AverageMass * RnaAverageO +
                PeriodicTable.GetElement("N").AverageMass * RnaAverageN;
        }

        internal static ChemicalFormula GetFormula(bool useAveratideModel)
        {
            return useAveratideModel ? GetAveratideFormula(AveratideMass) : GetAveragineFormula(AveragineMass);
        }

        internal static ChemicalFormula GetAveragineFormula(double mass)
        {
            double averagines = mass / (double)AveragineMass;
            ChemicalFormula formula = new ChemicalFormula();
            formula.Add("C", (int)Math.Round(averagines * AverageC, 0));
            formula.Add("H", (int)Math.Round(averagines * AverageH, 0));
            formula.Add("O", (int)Math.Round(averagines * AverageO, 0));
            formula.Add("N", (int)Math.Round(averagines * AverageN, 0));
            formula.Add("S", (int)Math.Round(averagines * AverageS, 0));
            return formula;
        }

        internal static ChemicalFormula GetAveratideFormula(double mass)
        {
            double averatides = mass / (double)AveratideMass;
            ChemicalFormula formula = new ChemicalFormula();
            formula.Add("C", (int)Math.Round(averatides * RnaAverageC, 0));
            formula.Add("H", (int)Math.Round(averatides * RnaAverageH, 0));
            formula.Add("O", (int)Math.Round(averatides * RnaAverageO, 0));
            formula.Add("N", (int)Math.Round(averatides * RnaAverageN, 0));
            return formula;
        }
    }
}
