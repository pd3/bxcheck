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

#ifndef __BXHASH_H__
#define __BXHASH_H__

#include <stdio.h>
#include <stdint.h>
#include <htslib/khash.h>
#include "bxcheck.h"

KHASH_MAP_INIT_INT(int2int, int)
typedef khash_t(int2int) bxhash_t;

bxhash_t *bxhash_init(char *fname);
void bxhash_destroy(bxhash_t *hash);
int bxhash_size(bxhash_t *hash);

static inline uint8_t char2nt(char nt)
{
    if ( nt=='A' || nt=='a' ) return 0;
    if ( nt=='C' || nt=='c' ) return 1;
    if ( nt=='G' || nt=='g' ) return 2;
    if ( nt=='T' || nt=='t' ) return 3;
    return 4;
}
static inline uint8_t char2nt_c(char nt)
{
    if ( nt=='A' || nt=='a' ) return 3;
    if ( nt=='C' || nt=='c' ) return 2;
    if ( nt=='G' || nt=='g' ) return 1;
    if ( nt=='T' || nt=='t' ) return 0;
    return 4;
}
static inline uint32_t bxhash_seq2key(char *seq)
{
    uint32_t key = 0, i;
    for (i=0; i<16; i++)
    {
        int base = char2nt(seq[i]);
        if ( base>3 ) return -1;
        key |= base << (i*2);
    }
    return key;
}
static inline void bxhash_key2seq(uint32_t num, char *seq)
{
    int i;
    for (i=0; i<16; i++) seq[i] = "ACGT"[ (num >> (i*2)) & 0x3 ];
    seq[16] = 0;
}
static inline int bxhash_has_key(bxhash_t *hash, uint32_t key)
{
    khint_t k = kh_get(int2int, hash, key);
    return k == kh_end(hash) ? 0 : k + 1;
}
static inline int _bxhash_has_seq(bxhash_t *hash, char *seq, int complement, int reverse)
{
    uint32_t key = 0;
    int i;
    for (i=0; i<16; i++)
    {
        int base = reverse ? seq[-i] : seq[i];
        base = complement ? char2nt_c(base) : char2nt(base);
        if ( base>3 ) return 0;
        key |= base << (i*2);
    }
    return bxhash_has_key(hash, key);
}
#define bxhash_has_seq(hash,seq)    _bxhash_has_seq(hash,seq,0,0)
#define bxhash_has_seq_c(hash,seq)  _bxhash_has_seq(hash,seq,1,0)
#define bxhash_has_seq_r(hash,seq)  _bxhash_has_seq(hash,seq,0,1)
#define bxhash_has_seq_cr(hash,seq) _bxhash_has_seq(hash,seq,1,1)

static inline int bxhash_put_key(bxhash_t *hash, uint32_t key)
{
    int ret;
    khint_t k = kh_put(int2int, hash, key, &ret);
    if ( ret < 0 ) return -1;
    if ( ret==0 ) kh_val(hash, k)++;
    else kh_val(hash, k) = 1;
    return 0;
}
static inline int bxhash_put_seq(bxhash_t *hash, char *seq)
{
    uint32_t key = bxhash_seq2key(seq);
    if ( key < 0 ) return key;
    return bxhash_put_key(hash, key);
}

#endif

