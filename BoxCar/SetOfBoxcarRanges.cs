using System;
using System.Collections;


namespace BoxCar
{
    /// <summary>
    /// A SetOfBoxcarRanges is an ArrayList of BoxcarRange objects.
    /// The number of BoxcarRange objects stored in the SetOfBoxcarRanges depends on how many boxes were specified in the scan setup (ex. 12)
    /// In Program.cs, SetOfBoxcarRanges is used for storing information about the boxcar scan method.
    /// More than one SetOfBoxcarRanges is usually needed depending on the number of boxcar (msx) scans. (program was tested with 2 or 3)
    /// </summary>
    public class SetOfBoxcarRanges : IEnumerable
    {
        public ArrayList Set { get; set; }

        /// <summary>
        /// Constructor that initializes the ArrayList Set to the given arrayList input.
        /// </summary>
        /// <param name="arrayList"></param>
        public SetOfBoxcarRanges(ArrayList arrayList)
        {
            Set = arrayList;
        }

        /// <summary>
        /// Default constructor for the SetOfBoxcarRanges class.
        /// </summary>
        public SetOfBoxcarRanges()
        {
            Set = new ArrayList();
        }

        /// <summary>
        /// Add a BoxcarRange object at the end of the SetOfBoxcarRanges.
        /// </summary>
        /// <param name="boxcarRange"></param>
        public void AddBoxcarRange(BoxcarRange boxcarRange)
        {
            Set.Add(boxcarRange);
        }

        /// <summary>
        /// The number of objects in the SetOfBoxcarRanges
        /// </summary>
        /// <returns></returns> int number of objects
        public int Count()
        {
            return Set.Count;
        }

        /// <summary>
        /// Gets the BoxcarRange element at the specified index in the SetOfBoxcarRanges ArrayList.
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns> BoxcarRange at index
        public BoxcarRange ElementAt(int index)
        {
            return (BoxcarRange)Set[index];
        }

        /// <summary>
        /// Inserts the given BoxcarRange object at the specified index.
        /// Replaces the previous object at that index.
        /// </summary>
        /// <param name="boxcarRange"></param>
        /// <param name="index"></param>
        public void Insert(BoxcarRange boxcarRange, int index)
        {
            Set[index] = boxcarRange;
        }


        // methods for iterating over the set of boxcar ranges  

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (IEnumerator)GetEnumerator();
        }

        public SetOfBoxcarRangesEnum GetEnumerator()
        {
            return new SetOfBoxcarRangesEnum(Set);
        }

        public class SetOfBoxcarRangesEnum : IEnumerator
        {
            public ArrayList Set;
            int position = -1;

            public SetOfBoxcarRangesEnum(ArrayList list)
            {
                Set = list;
            }

            public bool MoveNext()
            {
                position++;
                return (position < Set.Count);
            }

            public void Reset()
            {
                position = -1;
            }

            object IEnumerator.Current
            {
                get
                {
                    return Current;
                }
            }

            public BoxcarRange Current
            {
                get
                {
                    try
                    {
                        return (BoxcarRange)Set[position]; 
                    }
                    catch (IndexOutOfRangeException)
                    {
                        throw new InvalidOperationException();
                    }
                }
            }
        }
    }
}
