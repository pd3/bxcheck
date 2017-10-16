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
#include <signal.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/khash.h>

typedef struct
{
    int exp, cnt;
}
cnt_t;

KHASH_MAP_INIT_STR(str2cnt, cnt_t)

typedef struct
{
    khash_t(str2cnt) *str2cnt;
    htsFile *in, *out;
    bam_hdr_t *bam_hdr;
    char *bx_fname, *in_fname, *out_fname;
}
args_t;

static void error(const char *format, ...)
{
    if ( !format )
    {
        printf("About: Extract reads from BAM file by BX. (Intended for debugging)\n");
        printf("Usage: bxreads [OPTIONS]\n");
        printf("Options:\n");
        printf("    -b, --barcodes <file.txt>       \n");
        printf("    -i, --input-bam <file.bam>      \n");
        printf("    -o, --output-bam <file.bam>     \n");
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

int keep_the_read(args_t *args, const bam1_t *rec)
{
    uint8_t *aux = bam_aux_get(rec, "BX");   
    if ( !aux ) return 0;
    const char *tag = bam_aux2Z(aux);
    if ( !tag ) return 0;

    if ( strlen(tag)<16 ) error("really?? [%s]\n", tag);
    char tmp = tag[16];
    ((char*)tag)[16] = 0;
    khint_t k = kh_get(str2cnt, args->str2cnt, tag);
    ((char*)tag)[16] = tmp;

    if ( k == kh_end(args->str2cnt) ) return 0;
    kh_val(args->str2cnt, k).cnt++;

    return 1;
}

void init(args_t *args)
{
    FILE *fp = fopen(args->bx_fname,"r");
    if ( !fp ) error("Failed to open for reading: %s\n", args->bx_fname);
    args->str2cnt = kh_init(str2cnt);
    int cnt, bx, i;
    char tag[17];
    while ( fscanf(fp, "%d\t%d\n", &cnt,&bx)==2 )
    {
        for (i=0; i<16; i++)
        {
            int j = (bx >> (i*2)) & 0x3;
            tag[i] = "ACGT"[j];
        }
        tag[16] = 0;
        int ret;
        khint_t k = kh_put(str2cnt, args->str2cnt, strdup(tag), &ret);
        kh_val(args->str2cnt, k).exp = cnt;
        kh_val(args->str2cnt, k).cnt = 0;
    }
    if ( fclose(fp)!=0 ) error("Failed to close: %s\n", args->bx_fname);

    if ( (args->in = sam_open(args->in_fname,"r"))==0 ) error("Failed to open for reading: %s\n", args->in_fname);
    if ( (args->out = sam_open(args->out_fname,"wb"))==0 ) error("Failed to open for writing: %s\n", args->out_fname);
    args->bam_hdr = sam_hdr_read(args->in);
    if ( sam_hdr_write(args->out, args->bam_hdr) < 0 ) error("Could not write the bam header\n");
}

void destroy(args_t *args)
{
    bam_hdr_destroy(args->bam_hdr);
    if ( hts_close(args->in)!=0 ) error("Close failed: %s\n",args->in_fname);
    if ( hts_close(args->out)!=0 ) error("Close failed: %s\n",args->out_fname);

    khint_t k;
    for (k = 0; k < kh_end(args->str2cnt); ++k)
    {
        if ( !kh_exist(args->str2cnt,k) ) continue;
        int cnt = kh_val(args->str2cnt,k).cnt;
        int exp = kh_val(args->str2cnt,k).exp;
        if ( exp!=cnt ) fprintf(stderr,"Warning: expected %d reads, found %d\n", exp,cnt);
        free((char*)kh_key(args->str2cnt,k));
    }
    kh_destroy(str2cnt, args->str2cnt);
    free(args);
}

int main(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"barcodes", required_argument, NULL, 'b'},
        {"input-bam", required_argument, NULL, 'i'},
        {"output-bam", required_argument, NULL, 'o'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ( (opt=getopt_long(argc,argv,"?hb:i:o:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case 'i': 
                args->in_fname = optarg;
                break;
            case 'o': 
                args->out_fname = optarg;
                break;
            case 'b': 
                args->bx_fname = optarg;
                break;
            case '?':
            case 'h': error(NULL);
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( !args->in_fname || !args->out_fname || !args->bx_fname ) error(NULL);

    init(args);

    bam1_t *bam_line = bam_init1();
    while (sam_read1(args->in, args->bam_hdr, bam_line) >= 0)
    {
        if ( !keep_the_read(args, bam_line) ) continue;
        if ( sam_write1(args->out, args->bam_hdr, bam_line) < 0 ) error("Writing failed\n");
    }
    bam_destroy1(bam_line);

    destroy(args);

    return 0;
}

