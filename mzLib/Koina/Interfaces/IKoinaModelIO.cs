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
        public int MaxBatchSize { get; }
        public List<Dictionary<string, object>> ToBatchedRequests();
        public Task RunInferenceAsync();
    }
}
