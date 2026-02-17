using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PredictionClients.Koina.Interfaces
{
    public interface IPredictor<TModelInput, TModelOutput>
    {
        public Task<List<TModelOutput>> PredictAsync(List<TModelInput> inputs);
    }
}
