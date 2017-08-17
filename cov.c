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
#include <assert.h>
#include <htslib/hts.h>
#include "cov.h"
#include "rbuf.h"
#include "dist.h"

struct _cov_t
{
    rbuf_t rbuf;    // metadata for the round buffer pileup
    uint32_t *plp;  // pileup data
    int
        rid,        // previous rid, -1 on init
        off;        // first idx in rbuf corresponds to this coordinate (0-based)
    dist_t *dist;   // coverage histogram
};

cov_t *cov_init()
{
    cov_t *cov = (cov_t*) calloc(1,sizeof(cov_t));
    cov->dist = dist_init(4);
    cov->rid  = -1;
    return cov;
}

void cov_destroy(cov_t *cov)
{
    free(cov->plp);
    dist_destroy(cov->dist);
    free(cov);
}

void _cov_flush(cov_t *cov)
{
    int iplp;
    while ( (iplp=rbuf_shift(&cov->rbuf)) >=0 )
    {
        if ( cov->plp[iplp] ) dist_insert(cov->dist, cov->plp[iplp]);
        cov->plp[iplp] = 0;
    }
    cov->off = 0;
}

void cov_insert(cov_t *cov, int32_t rid, int32_t beg, int32_t end)
{
    assert( beg <= end );

    int plp_size = (end - beg + 1)*10;
    if ( cov->rbuf.m < plp_size )
    {
        if ( !cov->rbuf.m )
        {
            rbuf_init(&cov->rbuf, plp_size);
            cov->plp = (uint32_t*) calloc(plp_size, sizeof(*cov->plp));
        }
        else
            rbuf_expand0(&cov->rbuf, uint32_t, plp_size, cov->plp);
    }

    if ( rid != cov->rid )
    {
        _cov_flush(cov);
        cov->rid = rid;
    }

    int iplp = -1, ioff = 0;
    while ( rbuf_next(&cov->rbuf,&iplp) )
    {
        assert( cov->off + ioff <= beg );   // unsorted!
        if ( cov->off + ioff == beg ) break;
        if ( cov->plp[iplp] ) dist_insert(cov->dist, cov->plp[iplp]);
        cov->plp[iplp] = 0;
        ioff++;
    }
    if ( ioff ) 
    {
        rbuf_shift_n(&cov->rbuf, ioff);
        cov->off += ioff;
    }

    // debugging sanity checks
    if ( !cov->rbuf.n ) 
    {
        cov->off = beg;
        assert( cov->plp[cov->rbuf.f]==0 );
    }
    else
        assert( cov->off <= beg );

    int ibeg, iend, i;
    if ( end - cov->off >= cov->rbuf.n ) rbuf_append_n(&cov->rbuf, end - cov->off + 1);
    ibeg = rbuf_kth(&cov->rbuf, beg - cov->off);
    iend = rbuf_kth(&cov->rbuf, end - cov->off);
    if ( ibeg > iend )
    {
        for (i=ibeg; i<cov->rbuf.m; i++) cov->plp[i]++;
        ibeg = 0;
    }
    for (i=ibeg; i<=iend; i++) cov->plp[i]++;
}

int cov_n(cov_t *cov)
{
    if (cov->rbuf.n ) _cov_flush(cov);
    return dist_n(cov->dist);
}

uint64_t cov_get(cov_t *cov, uint32_t idx, uint32_t *beg, uint32_t *end)
{
    if ( cov->rbuf.n ) _cov_flush(cov);
    return dist_get(cov->dist, idx, beg, end);
}

