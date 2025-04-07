using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    public interface ISerializableIndexer
    {
        public void ClearIndex();
        public void SerializeIndex();
        public void DeserializeIndex();
    }
}
