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

#ifndef __COV_H__
#define __COV_H__

#include <stdio.h>
#include <inttypes.h>

typedef struct _cov_t cov_t;

cov_t *cov_init(int bin_size);
void cov_destroy(cov_t *cov);

void cov_insert(cov_t *cov, int32_t rid, int32_t beg, int32_t end);

int cov_nrid(cov_t *cov);                   // number of chromosomes
int cov_nbin(cov_t *cov, int32_t rid);      // number of bins on the chromosome

/*
   cov_get() 
   @idx:        from the interval [0,cov_n-1]
   @beg,end:    [beg,end)
 */
uint16_t cov_get(cov_t *cov, int32_t rid, int32_t bin, uint32_t *beg, uint32_t *end);

#endif

