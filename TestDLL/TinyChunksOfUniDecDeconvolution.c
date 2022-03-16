// Goal is to create tiny little chunks of UniDec's algorithm so that each component is completely tested and functional

#include "TinyChunksOfUniDecDeconvolution.h"

int MemoryAllocationOfBarr(void) {
    char* barr = NULL;

    int mlength, zlength, numclose,
        * mind = NULL,
        * zind = NULL,
        * closemind = NULL,
        * closezind = NULL,
        * closeind = NULL,
        * starttab = NULL,
        * endtab = NULL;
    float
        * mdist = NULL,
        * zdist = NULL,
        * mzdist = NULL,
        * rmzdist = NULL,
        * oldblur = NULL,
        * closeval = NULL,
        * closearray = NULL;

    free(mzdist);
    free(rmzdist);
    free(closeval);
    free(closearray);
    free(closemind);
    free(closezind);
    free(endtab);
    free(starttab);

    free(mdist);
    free(mind);
    free(zind);
    free(zdist);
    free(barr);
    free(closeind);
    return 1; 
}
