using PredictionClients.Koina.AbstractClasses;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PredictionClients.Koina.Interfaces
{
    public interface IPredictor<TModelInput, TModelOutput>
    {
        #region Properties
        List<TModelInput> ModelInputs { get; }
        bool[] ValidInputsMask { get; }
        List<TModelOutput> Predictions { get; }
        #endregion

        #region Methods
        List<TModelOutput> Predict(List<TModelInput> inputs);
        #endregion
    }
}
