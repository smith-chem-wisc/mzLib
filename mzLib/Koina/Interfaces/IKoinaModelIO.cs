using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Koina.Interfaces
{
    internal interface IKoinaModelIO
    {
        public string ModelName { get; }
        public int BatchSize { get; }
        public Dictionary<string, object> ToRequest(string? requestID = null);
        public void RunInference();
    }
}
