using Easy.Common.Extensions;
using SharpLearning.InputOutput.Csv;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;




namespace FlashLFQ.Alex_project
{
    public class ParallelSearch
    {

        readonly int threadNumber = Environment.ProcessorCount-1;
        readonly int[] threads;
        Dictionary<int, XIC> xicDict;
        XICGroups[] xicGroups;

        public ParallelSearch(Dictionary<int, XIC> xicDict)
        {
            threads = Enumerable.Range(0, threadNumber).ToArray();
            this.xicDict = xicDict;
            xicGroups = new XICGroups[threadNumber];
        }

        //public void run()
        //{
        //    Parallel.ForEach(threads, (currentThread) =>
        //    {
        //        List<XIC> xics = GroupedXIC(xicDict, currentThread);
        //        if (xics.Count() != 0)
        //        {
        //            xicGroups[currentThread] = new XICGroups(xics, 0.5, 0.1);
        //            draw(currentThread, xicGroups[currentThread]);
        //        }              
                
        //    });

        //    Console.WriteLine("The total sum is ");

        //}

        /// <summary>
        /// Try to group the XICs from a Big XIC data set. If we met the reference XIC, we will build the XICGroups.
        /// </summary>
        /// <param name="xicDict"></param>
        public static List<XIC> GroupedXIC(Dictionary<int, XIC> xicDict, int thread)
        {
            Dictionary<int, XIC> xicsToGroup = new Dictionary<int, XIC>();

            lock (xicDict) 
            {
                foreach (var xic in xicDict)
                {
                    if (xic.Value.Reference == true && xicsToGroup.Where(p => p.Value.Reference).Count() > 0)
                    {
                        break;
                    }

                    else
                    {
                        xicsToGroup.Add(xic.Key, xic.Value);
                        xicDict.Remove(xic.Key);
                    }
                }
            }

            return xicsToGroup.Select(p => p.Value).ToList();

        }


    }
}
