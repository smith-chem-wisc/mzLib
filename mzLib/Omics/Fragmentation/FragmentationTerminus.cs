﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.Fragmentation
{
    public enum FragmentationTerminus
        {
            Both, //N- and C-terminus
            N, //N-terminus only
            C, //C-terminus only
            None //used for internal fragments, could be used for top down intact mass?
        }
    
}
