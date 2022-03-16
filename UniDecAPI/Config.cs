using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
    [StructLayout(LayoutKind.Sequential)]
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

        public int imflag;
        //IM Parameters
        public float dtsig;
        public float csig;
        public float ccsub;
        public float ccslb;
        public float ccsbins;
        public float temp;
        public float press;
        public float volt;
        public float tcal1;
        public float tcal2;
        public float tcal3;
        public float tcal4;
        public float twaveflag;
        public float hmass;
        public float to;
        public float len;
        public float edc;
        public float nativeccsub;
        public float nativeccslb;
        public int baselineflag;
        public int noiseflag;
        public int zout;
        public int metamode;
        public float minmz;
        public float maxmz;
        public int mzbins;
        public float bsub;
        public float datareduction;
        public float peakwin;
        public float peakthresh;
        public float exwindow;
        public int exchoice;
        public int exchoicez;
        public float exthresh;
        public int exnorm;
        public int exnormz;
        public int peaknorm;
        public int orbimode;
        public int datanorm;
        //Experimental Parameters
        public int filterwidth;
        public float zerolog;
        public int lengthmz; // probably can get rid of 
        //public int mfilelen;
        public int isolength;
        // DoubleDec Parameters
        public int doubledec;
    }
}
