using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IsobaricQuant;
using NUnit.Framework;

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
    }
}

