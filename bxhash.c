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

#include <errno.h>
#include "bxhash.h"

int bxhash_size(bxhash_t *hash)
{
    return hash ? kh_size(hash) : 0; 
}

bxhash_t *bxhash_init(char *fname)
{
    bxhash_t *hash = kh_init(int2int);
    if ( !fname ) return hash;

    fprintf(stderr,"Processing known barcode sequences...\n");

    FILE *fp = fopen(fname,"r");
    if ( !fp ) error("%s: %s\n", fname,strerror(errno));
    char buf[18];
    int n = 0, i;
    while ( fread(buf,17,1,fp)==1 )
    {
        buf[17] = 0;
        n++;
        if ( buf[16]!='\n' ) error("Could not parse %d-th line in %s: %s\n", n, fname, buf);
        uint32_t bx = 0;
        for (i=0; i<16; i++)
        {
            int base = char2nt(buf[i]);
            if ( base>3 ) error("Could not parse %d-th line in %s: %s\n", n, fname, buf);
            bx |= base << (i*2);
        }
        bxhash_put_key(hash, bx);
    }
    if ( ferror(fp)!=0 ) error("Error reading %s\n", fname);
    if ( fclose(fp)!=0 ) error("close failed: %s\n", fname);
    fprintf(stderr,"Processed %d barcode sequences from %s\n", kh_size(hash), fname);
    return hash;
}

void bxhash_destroy(bxhash_t *hash)
{
    if ( hash ) kh_destroy(int2int, hash);
}

