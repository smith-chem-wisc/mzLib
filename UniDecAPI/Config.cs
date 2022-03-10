using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public struct Config
	{
        public int numit;// number iterations
        public int numz; // number of charge states
        public int endz; // final charge state
        public int startz; // initial charge state
        public float zsig; //
        public float psig;
        public float beta;
        public float mzsig;
        public float msig;
        public float molig;
        public float massub;
        public float masslb;
        public int psfun;
        public float mtabsig;
        public int mflag;
        public float massbins;
        public int limitflag;
        public float psthresh;
        public int speedyflag;
        public int linflag;
        public int aggressiveflag;
        public float adductmass;
        public int rawflag;
        public float nativezub;
        public float nativezlb;
        public int poolflag;
        public int manualflag;
        public float intthresh;
        public float peakshapeinflate;
        public float killmass;
        public int fixedmassaxis;
        public int isotopemode;
        public int filetype;
        public int imflag;
    }
}
