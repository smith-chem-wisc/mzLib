using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public interface IMs1FeatureFile
    {
        public IEnumerable<ISingleChargeMs1Feature> GetMs1Features();
    }
}
