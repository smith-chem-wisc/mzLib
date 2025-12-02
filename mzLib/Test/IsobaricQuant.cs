using IsobaricQuant;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    [TestFixture]
    internal class IsobaricQuant
    {
        [Test]
        public void Ctor_WithValidInput()
        {
            var input = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>
            {
                (2, new[]{49,56,30}),
                (3, new[]{86,45,97}),
                (6, new[]{61,59,44}),
                (1, new[]{35,32,86}),
                (2, new[]{35,31,13}),
                (6, new[]{49,74,18}),
                (4, new[]{87,9,23}),
                (4, new[]{12,63,63}),
                (5, new[]{0,84,96}),
                (3, new[]{54,53,54}),
                (1, new[]{24,11,9}),
                (5, new[]{35,25,77})
            };
            var inputList = new List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> { input };
            IsobaricQuantEngine ibq = new IsobaricQuantEngine(inputList);
            var o1 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);
            var o2 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.Max);
            var o3 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.Median);
            var o4 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN, topN: 2);
            var o5 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN, lowFraction: 0.2);
            var o6 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN, normalize: false);
            var o7 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN, topN: 2, referenceChannel: 3);
            var input2 = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>
            {
                (7, new[]{30,56,49}),
                (3, new[]{97,45,86}),
                (6, new[]{44,59,61}),
                (1, new[]{86,32,35}),
                (8, new[]{13,31,35}),
                (6, new[]{18,74,49}),
                (4, new[]{23,9,87}),
                (4, new[]{63,63,12}),
                (9, new[]{96,84,0}),
                (3, new[]{54,53,54}),
                (1, new[]{9,11,24}),
                (5, new[]{77,25,35})
            };
            inputList.Add(input2);
            IsobaricQuantEngine ibq2 = new IsobaricQuantEngine(inputList);
            var o8 = ibq2.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);
        }
        [Test]
        public void ProcessPsms()
        {
            string psmFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\2025-12-02-13-52-34\Task1-SearchTask\AllPSMs.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            var input = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>();
            int counter = 0;
            foreach (var psm in parsedPsms)
            {
                input.Add((peptideFullSequenceHash: psm.FullSequence.GetHashCode(), reporterIntensities: RandomIntegerArray(6, 50, 100)));
                counter++;
                if (counter >= 20000)
                {
                    break;
                }
            }
            var inputList = new List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> { input };
            IsobaricQuantEngine ibq = new IsobaricQuantEngine(inputList);
            var o1 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);
            
            psmFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\2025-12-02-13-54-13\Task1-SearchTask\AllPSMs.psmtsv";
            parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out warnings);
            input = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>();
            counter = 0;
            foreach (var psm in parsedPsms)
            {
                input.Add((peptideFullSequenceHash: psm.FullSequence.GetHashCode(), reporterIntensities: RandomIntegerArray(6, 50, 100)));
                counter++;
                if (counter >= 20000)
                {
                    break;
                }
            }
            inputList = new List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> { input };
            var ibq2 = new IsobaricQuantEngine(inputList);
            var o2 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);
            int[] newChannelLabels1 = new int[]
            {
                "A549_1".GetHashCode(),
                "A549_2".GetHashCode(),
                "A549_3".GetHashCode(),
                "A549_4".GetHashCode(),
                "A549_5".GetHashCode(),
                "A549_6".GetHashCode()
            };
            var r1 = ibq.ReplaceChannelIndexBySampleKey(o1, newChannelLabels1);
            var r2 = ibq2.ReplaceChannelIndexBySampleKey(o2, newChannelLabels1);
            var g = ibq.GlobalResults(new List<ConcurrentDictionary<int, ConcurrentDictionary<int, int>>> { r1, r2 });
        }
        [Test]
        public void ProcessProteinPsms()
        {
            string psmFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\2025-12-02-13-52-34\Task1-SearchTask\AllPSMs.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            var input = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>();
            int counter = 0;
            foreach (var psm in parsedPsms)
            {
                input.Add((peptideFullSequenceHash: psm.ProteinAccession.GetHashCode(), reporterIntensities: RandomIntegerArray(6, 50, 100)));
                counter++;
                if (counter >= 20000)
                {
                    break;
                }
            }
            var inputList = new List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> { input };
            IsobaricQuantEngine ibq = new IsobaricQuantEngine(inputList);
            var o1 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);

            psmFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\2025-12-02-13-54-13\Task1-SearchTask\AllPSMs.psmtsv";
            parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out warnings);
            input = new List<(int peptideFullSequenceHash, int[] reporterIntensities)>();
            counter = 0;
            foreach (var psm in parsedPsms)
            {
                input.Add((peptideFullSequenceHash: psm.ProteinAccession.GetHashCode(), reporterIntensities: RandomIntegerArray(6, 50, 100)));
                counter++;
                if (counter >= 20000)
                {
                    break;
                }
            }
            inputList = new List<List<(int peptideFullSequenceHash, int[] reporterIntensities)>> { input };
            var ibq2 = new IsobaricQuantEngine(inputList);
            var o2 = ibq.Process(aggregateType: IsobaricQuantEngine.AggregateType.SumTopN);
            int[] newChannelLabels1 = new int[]
            {
                "A549_1".GetHashCode(),
                "A549_2".GetHashCode(),
                "A549_3".GetHashCode(),
                "A549_4".GetHashCode(),
                "A549_5".GetHashCode(),
                "A549_6".GetHashCode()
            };
            var r1 = ibq.ReplaceChannelIndexBySampleKey(o1, newChannelLabels1);
            var r2 = ibq2.ReplaceChannelIndexBySampleKey(o2, newChannelLabels1);
            var g = ibq.GlobalResults(new List<ConcurrentDictionary<int, ConcurrentDictionary<int, int>>> { r1, r2 });
        }
        public int[] RandomIntegerArray(int length, int minValue, int maxValue)
        {
            Random rand = new Random();
            int[] array = new int[length];
            for (int i = 0; i < length; i++)
            {
                array[i] = rand.Next(minValue, maxValue);
            }
            return array;
        }
    }
}

