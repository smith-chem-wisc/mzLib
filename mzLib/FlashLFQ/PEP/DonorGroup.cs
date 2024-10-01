using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.PEP
{
    public class DonorGroup : IEnumerable<ChromatographicPeak>
    {
        public Identification DonorId { get; }
        public List<ChromatographicPeak> TargetAcceptors { get; }
        public List<ChromatographicPeak> DecoyAcceptors { get; }

        public DonorGroup(Identification donorId, List<ChromatographicPeak> targetAcceptors, List<ChromatographicPeak> decoyAcceptors)
        {
            DonorId = donorId;
            TargetAcceptors = targetAcceptors;
            DecoyAcceptors = decoyAcceptors;
        }

        public double BestTargetMbrScore => TargetAcceptors.Count == 0 ? 0 : TargetAcceptors.Max(acceptor => acceptor.MbrScore);

        public IEnumerator<ChromatographicPeak> GetEnumerator()
        {
            return TargetAcceptors.Concat(DecoyAcceptors).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

    }
}
