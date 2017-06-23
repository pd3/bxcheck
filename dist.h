/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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

#ifndef __DIST_H__
#define __DIST_H__

#include <stdio.h>
#include <inttypes.h>

typedef struct _dist_t dist_t;

/*
 *  dist_init() - init bins
 */
dist_t *dist_init(int npow);
void dist_destroy(dist_t *dist);

/*
    dist_n() - get the number of bins
 */
int dist_n(dist_t *dist);

/*
    dist_insert() - insert new value
 */
uint32_t dist_insert(dist_t *dist, uint32_t value);

/*
   dist_get() 
   @idx:        from the interval [0,dist_n-1]
   @beg,end:    [beg,end)
 */
uint64_t dist_get(dist_t *dist, uint32_t idx, uint32_t *beg, uint32_t *end);

#endif
