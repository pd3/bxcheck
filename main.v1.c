/*  bxcheck.c -- Collect chromium BX stats

    Copyright (C) 2017 Genome Research Ltd.

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  
*/

#include <unistd.h>     // for isatty()
#include <inttypes.h>   // for PRIu64
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include "kheap.h"

/*
    - what is the number of barcodes?
    - per barcode, what is distribution of number of read pairs?
    - what is the distribution of estimated fragment size?
    - what is the distribution of number of fragments per barcode?
    - what is the typical spacing of read pairs within a fragment?
    - what is average fragment depth; how many fragments covering each position in the genome (distribution?)?
*/

/*
  For counting the number of barcodes, lives throughout the whole bam.
  The values are 1=active, 2=done
*/
KHASH_MAP_INIT_STR(str2uint8, uint8_t)


/* 
  For estimating the number of fragments, size distribution and spacing within a barcode.
  Short-lived, controlled by --max-fragment-size, contains the active str2bx barcodes
*/
typedef struct
{
    uint32_t npos, mpos, *pos;
    const char *str;    // the BX tag
}
bx_t;
KHASH_MAP_INIT_STR(str2bx, bx_t*)

static inline int bx_is_smaller(bx_t **a, bx_t **b)
{
    return (*a)->pos[0] < (*b)->pos[0] ? 1 : 0; 
}
KHEAP_INIT(bx, bx_t*, bx_is_smaller)
typedef khp_bx_t bx_heap_t;


static inline unsigned hash_Wang(unsigned key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

static inline unsigned hash_X31_Wang(const char *s)
{
    unsigned h = *s;
    if (h) {
        for (++s ; *s; ++s) h = (h << 5) - h + *s;
        return hash_Wang(h);
    } else return 0;
}

typedef struct {
    unsigned key;
    uint32_t pos.....
} elem_t;

static inline int elem_lt(elem_t x, elem_t y)
{
    if (x.key < y.key) return 1;
    if (x.key == y.key) {
        int t;
        t = strcmp(bam_get_qname(x.b), bam_get_qname(y.b));
        if (t < 0) return 1;
        return (t == 0 && ((x.b->core.flag>>6&3) < (y.b->core.flag>>6&3)));
    } else return 0;
}

KSORT_INIT(bamshuf, elem_t, elem_lt)



typedef struct
{
    uint32_t tid, pos;
    uint64_t nanomalous;            // the number of reads outside of --max-fragment-size
    uint64_t nreads, nprocessed;    // the number reads total and the number of reads with BX tag
    int max_fragment_size;          // maximum barcode length, reads falling outside are counted in nanomalous
    khash_t(str2uint8) *bx_all;     // unique barcodes and their status (1=active, 2=done), lives throughout the whole bam
    khash_t(str2bx) *bx_active;     // hash of active barcodes
    bx_heap_t *bx_heap;             // heap of active barcodes
    htsFile *bam_fh;
    bam_hdr_t *bam_hdr;
    int argc;
    char **argv;
}
args_t;

static void error(const char *format, ...)
{
    if ( !format )
    {
        printf("About: The program collects BX statistics from BAM files\n");
        printf("Usage: bxcheck [OPTIONS] file.bam\n");
        printf("Options:\n");
        printf("    -M, --max-fragment-size <float>      Maximum gap between two fragments [10e6]\n");
        printf("\n");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(-1);
}

void flush(args_t *args, uint32_t pos)
{
    bx_heap_t *heap = args->bx_heap;
    while ( heap->ndat )
    {
        bx_t *bx = heap->dat[0];
        if ( pos!=UINT32_MAX && pos - bx->pos[0] < args->max_fragment_size ) return;

        printf("BX\t%s\t%u\t%u", bx->str, bx->npos, bx->pos[0]);
        int i;
        for (i=1; i<bx->npos; i++) printf("\t%u", bx->pos[i] - bx->pos[i-1]);
        printf("\n");

        khint_t k = kh_get(str2bx, args->bx_active, bx->str);
        kh_del(str2bx, args->bx_active, k);

        k = kh_get(str2uint8, args->bx_all, bx->str);
        kh_val(args->bx_all, k) = 2;         // mark this barcode as done

        free(bx->pos);
        free(bx);
        khp_delete(bx, heap);
    }
}

void process(args_t *args, const bam1_t *rec)
{
    args->nreads++;
    uint8_t *aux = bam_aux_get(rec,"BX");
    if ( !aux ) return;
    const char *tag = bam_aux2Z(aux);
    if ( !tag ) return;

    if ( args->nprocessed )
        flush(args, args->tid==rec->core.tid ? rec->core.pos : UINT32_MAX);

    args->nprocessed++;
    args->tid = rec->core.tid;
    args->pos = rec->core.pos;

    khint_t k = kh_get(str2uint8, args->bx_all, tag);
    if ( k == kh_end(args->bx_all) )
    {
        int ret;
        tag = strdup(tag);
        k = kh_put(str2uint8, args->bx_all, tag, &ret);
        if ( ret < 0 ) error("Failed to store the barcode: %s\n", tag);
        kh_val(args->bx_all, k) = 0;
    }
    else
        tag = kh_key(args->bx_all, k);   // get a stable, long-lived copy of the string

    int status = kh_val(args->bx_all, k);
    if ( status == 2 )
    {
        args->nanomalous++;
        return;
    }

    bx_t *bx = NULL;
    if ( status == 0 )
    {
        int ret;
        kh_val(args->bx_all, k) = 1;
        bx = (bx_t*) calloc(1,sizeof(bx_t));
        bx->str = tag;
        k = kh_put(str2bx, args->bx_active, tag, &ret);
        kh_val(args->bx_active, k) = bx;
    }
    else
    {
        k = kh_get(str2bx, args->bx_active, tag);
        assert( k != kh_end(args->bx_active) );
        bx = kh_val(args->bx_active, k);
    }

    bx->npos++;
    hts_expand(uint32_t, bx->npos, bx->mpos, bx->pos);
    bx->pos[bx->npos-1] = rec->core.pos;

    assert( bx->npos==1 || rec->core.pos - bx->pos[bx->npos-2] < args->max_fragment_size );

    if ( status == 0 )
        khp_insert(bx, args->bx_heap, &bx);
}

void output(args_t *args)
{
    printf("SN\tn_barcodes\t%d\n", kh_size(args->bx_all));
    printf("SN\tn_reads\t%"PRIu64"\n", args->nreads);
    printf("SN\tn_processed\t%"PRIu64"\n", args->nprocessed);
    printf("SN\tn_anomalous\t%"PRIu64"\n", args->nanomalous);
}

void init(args_t *args)
{
    args->bx_all = kh_init(str2uint8);
    args->bx_active = kh_init(str2bx);
    args->bx_heap = khp_init(bx);
}

void destroy(args_t *args)
{
    khint_t k;
    for (k = 0; k < kh_end(args->bx_all); k++)
        if ( kh_exist(args->bx_all, k) ) free((char*)kh_key(args->bx_all, k));
    khp_destroy(bx, args->bx_heap);
    kh_destroy(str2bx, args->bx_active);
    kh_destroy(str2uint8, args->bx_all);
    free(args);
}

int main(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->max_fragment_size = 10e6;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"max-fragment-gap", required_argument, NULL, 1},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    int opt;
    while ( (opt=getopt_long(argc,argv,"?hM:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case 'f': 
                args->max_fragment_size = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --max-fragment-size %s\n", optarg);
                break;
            case '?':
            case 'h': error(NULL);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *bam_fname = NULL;
    if ( optind<argc )
        bam_fname = argv[optind++];
    if ( !bam_fname )
    {
        if ( isatty(STDIN_FILENO) ) error(NULL);
        bam_fname = "-";
    }

    if ((args->bam_fh = sam_open(bam_fname, "r")) == 0) error("Failed to open: %s\n", bam_fname);
    args->bam_hdr = sam_hdr_read(args->bam_fh);
    bam1_t *bam_line = bam_init1();

    init(args);
fprintf(stderr,"[todo]: filter reads by flag\n");

    while (sam_read1(args->bam_fh, args->bam_hdr, bam_line) >= 0) process(args, bam_line);
    flush(args, UINT32_MAX);

    output(args);

    bam_destroy1(bam_line);
    bam_hdr_destroy(args->bam_hdr);
    if ( hts_close(args->bam_fh)!=0 ) error("Close failed: %s\n",strcmp("-",bam_fname)?bam_fname:"STDIN");
    destroy(args);

    return 0;
}
