/* The MIT License

   Copyright (c) 2017 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>
   
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include "cov.h"

typedef struct
{
    int nbin, mbin;
    uint16_t *bin;
}
bins_t;

struct _cov_t
{
    int bin_size;
    int nrid;
    bins_t *dat;
};

cov_t *cov_init(int bin_size)
{
    cov_t *cov = (cov_t*) calloc(1,sizeof(cov_t));
    cov->bin_size = bin_size;
    return cov;
}

void cov_destroy(cov_t *cov)
{
    int i;
    for (i=0; i<cov->nrid; i++) free(cov->dat[i].bin);
    free(cov->dat);
    free(cov);
}

void cov_insert(cov_t *cov, int32_t rid, int32_t beg, int32_t end)
{
    int i;
    if ( rid+1 > cov->nrid )
    {
        cov->dat = (bins_t*) realloc(cov->dat, sizeof(bins_t)*(rid+1));
        for (i=cov->nrid; i<=rid; i++)
            memset(cov->dat + cov->nrid, 0, sizeof(bins_t)*(rid + 1 - cov->nrid));
        cov->nrid = rid+1;
    }
    bins_t *dat = &cov->dat[rid];
    int ibeg = beg / cov->bin_size;
    int iend = end / cov->bin_size;
    hts_expand0(uint16_t, iend+1, dat->mbin, dat->bin);
    for (i=ibeg; i<=iend; i++) 
        if ( dat->bin[i] < 65535 ) dat->bin[i]++;
    if ( iend >= dat->nbin ) dat->nbin = iend+1;
}

int cov_nrid(cov_t *cov) { return cov->nrid; }
int cov_nbin(cov_t *cov, int32_t rid) { return cov->dat[rid].nbin; }

uint16_t cov_get(cov_t *cov, int32_t rid, int32_t bin, uint32_t *beg, uint32_t *end)
{
    if ( beg ) *beg = cov->bin_size * bin;
    if ( end ) *end = cov->bin_size * (bin + 1);
    return cov->dat[rid].bin[bin];
}


