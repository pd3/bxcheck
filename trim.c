/*  fq.c -- Extracts barcodes

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
/*
   The error correction algorithm:
    - build the list of known barcode sequences (the hash bx_known)
    - build the list of known barcode sequences found in the fq (the hash bx_fq)
    - take reads one by one and calculate the most likely barcode:
        P(BX|seq) = P(seq|BX) * P(BX) / sum_{bx} P(seq|bx) P(bx)
      where
        P(seq|BX) = \prod_i Q_i
        P(BX) = n(BX)/\sum_{bx} n(bx)
      and
        Q_i .. the Illumina quality of i-th base in the barcode sequence
        n(BX) .. the number of reads with the barcode BX
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include "bxcheck.h"
#include "bxhash.h"
#include "dist.h"

typedef struct
{
    uint64_t n, mkey, mqual;
    uint32_t *key;  // barcode keys, created with bxhash_seq2key, N's replaced with A's
    uint8_t *qual;  // see base_qual_set and base_qual_get, each BQ in 3bits, reduced by Illumina's binning 
}
reads_t;

typedef struct
{
    reads_t reads;
    bxhash_t        // multiplicity of reads in barcodes:
        *bx_known,  // .. the pool of all known barcodes
        *bx_fq,     // .. all barcodes encountered in the fastq, except if they contain N's
        *bx_out_kn, // .. barcodes after correction, known barcodes
        *bx_out_un; // .. barcodes after correction, unknown barcodes
    dist_t
        *rpb_raw_kn,    // number of reads in known raw barcodes
        *rpb_raw_un,    // number of reads in unknown raw barcodes
        *rpb_mod_ori,   // original number of reads in modifed barcodes
        *rpb_mod_new;   // new number of reads in modified barcodes
    uint64_t nfew_reads;    // reads that could be changed but were not because the new barcode had too few reads
    int argc, trust_nreads, good_bx_nreads, discard;
    char **argv, *prefix, *fq1, *fq2;
    kstring_t str;
    FILE *fh_log;
}
args_t;


static double qbin2prob[8];  // see base_qual_get16()
static double qbin2probc[8]; // see base_qual_get16()
static void base_qual_init_qbin2prob(void)
{
    // prob of the base being incorrect
    qbin2prob[0] = 0.75;             // N
    qbin2prob[1] = pow(10,-0.1* 6);  // Q = 2-9
    qbin2prob[2] = pow(10,-0.1*15);  //     10-19
    qbin2prob[3] = pow(10,-0.1*22);  //     20-24
    qbin2prob[4] = pow(10,-0.1*27);  //     25-29
    qbin2prob[5] = pow(10,-0.1*33);  //     30-34
    qbin2prob[6] = pow(10,-0.1*37);  //     35-39
    qbin2prob[7] = pow(10,-0.1*40);  //     >=40 

    // log prob of the complement - of a different, randomly selected base being correct
    int i;
    for (i=0; i<8; i++)
        qbin2probc[i] = log(qbin2prob[i]/3.);

    // prob of the base being correct
    for (i=0; i<8; i++)
        qbin2prob[i] = log(1 - qbin2prob[i]);
}

static void destroy(args_t *args)
{
    bxhash_destroy(args->bx_known);
    bxhash_destroy(args->bx_fq);
    fclose(args->fh_log);
    free(args->reads.key);
    free(args->reads.qual);
    free(args->str.s);
    free(args);
}

static inline void base_qual_set16(uint8_t *dst, uint8_t *quals)
{
    /*
        Quality score reduced into 3 bits, similar to Illumina's quality binning 
            http://www.uppmax.uu.se/digitalAssets/557/c_557912-l_1-k_whitepaper_datacompression.pdf

        The 16 barcode qualities are stored in 6 bytes, organised as follows:
            01234567|01234567|01234567|..   - bits
            00011122 23334445 55666777 ..   - i-th quality, 
    */
    const int qual2bin[] = {0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6};

    memset(dst, 0, 6);

    int ibase;
    for (ibase=0; ibase<16; ibase++)
    {
        int qbin = quals[ibase] >= 40 ? 7 : qual2bin[quals[ibase]];
        // which of the 6 bytes to use and how many - one or two?
        int idst = ibase * 3 / 8;
        int ibit = ibase * 3 % 8;
        if ( ibit < 6 )
            dst[idst] |= qbin << ibit;
        else if ( ibit==6 )
        {
            dst[idst] |= (qbin & 0x3) << ibit;  // low bits into the first byte
            dst[idst+1] |= (qbin & 0x4) >> 2;   // the high bit into the second byte
        }
        else    // ibit==7
        {
            dst[idst] |= (qbin & 0x1) << ibit;  // the low bit into the first byte
            dst[idst+1] |= (qbin & 0x6) >> 1;   // high bit into the second byte
        }
    }
}
static inline void base_qual_get16(double *prob, double *probc, uint8_t *qbin)
{
    int ibase;
    for (ibase=0; ibase<16; ibase++)
    {
        int ibin = ibase * 3 / 8;
        int ibit = ibase * 3 % 8;
        int bin;
        if ( ibit < 6 )
            bin = (qbin[ibin] >> ibit) & 0x7;
        else if ( ibit==6 )
        {
            bin = (qbin[ibin] >> ibit) & 0x3;
            bin |= (qbin[ibin] & 0x1) << 2;
        }
        else    // ibit==7
        {
            bin = (qbin[ibin] >> ibit) & 0x1;
            bin |= (qbin[ibin] & 0x3) << 1;
        }
        prob[ibase]  = qbin2prob[bin];
        probc[ibase] = qbin2probc[bin];
    }
}
static int add_read(reads_t *reads, char *seq, char *qual)
{
    reads->n++;
    hts_expand(uint32_t, reads->n, reads->mkey, reads->key);
    hts_expand(uint8_t, 6*reads->n, reads->mqual, reads->qual);

    uint8_t tmp[16], *quals = reads->qual + (reads->n-1)*6;
    int i, has_n = 0;
    for (i=0; i<16; i++)
    {
        if ( char2nt(seq[i]) > 3 ) tmp[i] = 0, seq[i] = 'A', has_n = 1;
        else tmp[i] = qual[i] - 33;
    }
    base_qual_set16(quals, tmp);

    reads->key[reads->n - 1] = bxhash_seq2key(seq);
    return has_n;
}

static void read_fq_bx(args_t *args)
{
    args->bx_fq = bxhash_init(NULL);

    kstring_t str = {0,0,0}, name = {0,0,0}, seq = {0,0,0}, qual = {0,0,0};
    htsFile *fp = hts_open(args->fq1,"r");
    if ( !fp ) error("Failed to read %s\n", args->fq1);
    int nskip = args->fq2 ? 8 : 4;
    while ( 1 )
    {
        //  @ST-E00143:240:HVTNFCCXX:2:1101:2514:1766 1:N:0:0
        //  TGAAAGATCTTGAAGCGAACAGATTAAGTCTGTTTCCAAATAAAAGGTT...
        //  +
        //  AAAFFJJJJJJJJJJJJJJJJFJFJJJFF<AJJJFJ-J<FJJJFJF<-F...

        if ( hts_getline(fp, KS_SEP_LINE, &name) <= 0 ) break;
        if ( name.s[0]!='@' ) error("Could not parse %s, line %d: %s\n", args->fq1, args->reads.n*nskip+1, name.s);
        if ( hts_getline(fp, KS_SEP_LINE, &seq) <= 0 ) error("Truncated file? %s\n", args->fq1);
        if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Truncated file? %s\n", args->fq1);
        if ( hts_getline(fp, KS_SEP_LINE, &qual) <= 0 ) error("Truncated file? %s\n", args->fq1);

        if ( !args->fq2 )   // interleaved fastq
        {
            if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Truncated file? %s\n", args->fq1);
            int i;
            for (i=0; i<name.l; i++) if ( isspace(name.s[i]) ) { name.s[i] = 0; break; }
            for (i=0; i<str.l; i++) if ( isspace(str.s[i]) ) { str.s[i] = 0; break; }
            if ( strcmp(name.s,str.s) ) error("Not an interleaved fastq file: %s vs %s, line %d\n",name.s,str.s,args->reads.n*8+1);
        }

        if ( add_read(&args->reads, seq.s, qual.s) == 0 )   // real BX sequence, does not contain N
            bxhash_put_key(args->bx_fq, args->reads.key[args->reads.n-1]);

        if ( !args->fq2 )   // interleaved fastq
        {
            if ( hts_getline(fp, KS_SEP_LINE, &seq) <= 0 ) error("Truncated file? %s\n", args->fq1);
            if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Truncated file? %s\n", args->fq1);
            if ( hts_getline(fp, KS_SEP_LINE, &qual) <= 0 ) error("Truncated file? %s\n", args->fq1);
        }
    }
    free(str.s);
    free(name.s);
    free(seq.s);
    free(qual.s);
    if ( hts_close(fp)!=0 ) error("Close failed: %s\n", args->fq1);
}

static inline double logsumexp2(double a, double b)
{
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}
static int correct_bx(args_t *args, int iread, uint32_t *fixed_key)
{
    uint32_t
        keys[64],       // single error-sequences observed in the data
        nbx[64],        // bx frequencies for the 64 possible single-error sequences
        nbx_tot = 0;    // \sum_i nbx[i]
    double
        bq[16],         // log P(base_is_correct)
        bqc[16],        // log P(different_base_is_correct)
        bq_prod = 0,
        bx_lk[64],      // P(seq|BX) likelihoods
        bx_prob[64];    // P(seq|BX)*P(BX)

    int i, j, k, n = 0;

    const uint32_t raw_key = args->reads.key[iread];
    int is_known = bxhash_has_key(args->bx_known, raw_key);

    // this misclassifies reads with N's
    k = bxhash_has_key(args->bx_fq, raw_key);
    uint32_t nori = 0;
    if ( k ) 
        nori = kh_val(args->bx_fq, k-1);

    if ( is_known )
    {
        dist_insert(args->rpb_raw_kn, nori);
        if ( nori >= args->trust_nreads ) 
        {
            *fixed_key = raw_key;
            bxhash_put_key(args->bx_out_kn, raw_key);
            return 1;
        }
    }
    else
        dist_insert(args->rpb_raw_un, nori);

    base_qual_get16(bq, bqc, args->reads.qual + iread*6);
    for (i=0; i<16; i++) bq_prod += bq[i];

    if ( is_known )
    {
        keys[n]  = raw_key;
        nbx[n]   = nori;
        nbx_tot += nori;
        bx_lk[n] = bq_prod;
        n++;
    }
    for (i=0; i<16; i++)
        for (j=0; j<4; j++)
        {
            uint32_t key = raw_key;
            key &= ~(0x3 << i*2);
            key |= j << i*2;
            if ( key == raw_key ) continue;
            if ( !bxhash_has_key(args->bx_known, key) ) continue;
            k = bxhash_has_key(args->bx_fq, key);
            if ( !k ) continue;
            keys[n]  = key;
            nbx[n]   = kh_val(args->bx_fq,k-1);
            nbx_tot += nbx[n];
            bx_lk[n] = bq_prod - bq[i] + bqc[i];    // replace P(base_is_correct) with P(different_base_is_correct)
            n++;
        }

    if ( n==0 ) // failed to find a known barcode with a single-nt difference
    {
        bxhash_put_key(args->bx_out_un, raw_key);
        *fixed_key = raw_key;
        return -1;
    }

    int imax = 0;
    double max, bx_prob_tot;
    max = bx_prob_tot = bx_prob[0] = bx_lk[0] + log(nbx[0]); 
    for (i=1; i<n; i++)
    {
        bx_prob[i] = bx_lk[i] + log(nbx[i]); 
        bx_prob_tot = logsumexp2(bx_prob[i], bx_prob_tot);
        if ( max < bx_prob[i] )
        {
            max  = bx_prob[i];
            imax = i;
        }
    }
    if ( keys[imax] != raw_key && nbx[imax] < args->good_bx_nreads ) // failed to find a barcode with enough reads
    {
        *fixed_key = raw_key;
        bxhash_put_key(is_known ? args->bx_out_kn : args->bx_out_un, raw_key);
        args->nfew_reads++;
        if ( is_known ) return 0;
        return -1;
    }
    *fixed_key = keys[imax];
    bxhash_put_key(args->bx_out_kn, *fixed_key);

    if ( *fixed_key != raw_key )
    {
        dist_insert(args->rpb_mod_ori, nori);
        dist_insert(args->rpb_mod_new, nbx[imax]);
    }

    return 0;
}

static void analyze_fq_bx(args_t *args)
{
    kstring_t str = {0,0,0}, name = {0,0,0}, name2 = {0,0,0}, out_name = {0,0,0}, seq = {0,0,0}, qual = {0,0,0};
    htsFile *in1 = hts_open(args->fq1,"r"), *in2 = args->fq2 ? hts_open(args->fq2,"r") : NULL;
    if ( !in1 ) error("Failed to read %s\n", args->fq1);
    if ( args->fq2 && !in2 ) error("Failed to read %s\n", args->fq2);
    ksprintf(&str,"%s.1.fq.gz",args->prefix);
    htsFile *out1 = hts_open(str.s,"wz");
    if ( !out1 ) error("Failed to open for writing %s: %s\n", str.s);
    char *out1_fname = strdup(str.s);
    str.l = 0;
    ksprintf(&str,"%s.2.fq.gz",args->prefix);
    htsFile *out2 = hts_open(str.s,"wz");
    if ( !out2 ) error("Failed to open for writing %s: %s\n", str.s);
    char *out2_fname = strdup(str.s);
    int nskip = args->fq2 ? 8 : 4;

    args->rpb_mod_ori = dist_init(4);
    args->rpb_mod_new = dist_init(4);
    args->rpb_raw_kn = dist_init(4);
    args->rpb_raw_un = dist_init(4);
    args->bx_out_kn = bxhash_init(NULL);
    args->bx_out_un = bxhash_init(NULL);
    uint64_t nchanged = 0, nfailed = 0, ntrusted = 0, nwr = 0, nwr_bx = 0, nwr_bx_mod = 0;
    int i,j, ret;
    for (i=0; i<args->reads.n; i++)
    {
        uint32_t key;
        ret = correct_bx(args, i, &key);
        if ( ret == -1  ) nfailed++;
        else if ( ret == 1 ) ntrusted++;
        else if ( args->reads.key[i] != key ) nchanged++;

        //  @ST-E00143:240:HVTNFCCXX:2:1101:2514:1766 1:N:0:0
        //  TGAAAGATCTTGAAGCGAACAGATTAAGTCTGTTTCCAAATAAAAGGTT...
        //  +
        //  AAAFFJJJJJJJJJJJJJJJJFJFJJJFF<AJJJFJ-J<FJJJFJF<-F...

        if ( hts_getline(in1, KS_SEP_LINE, &name) <= 0 ) break;
        if ( name.s[0]!='@' ) error("Could not parse %s, line %d: %s\n", args->fq1, args->reads.n*nskip+1, name.s);
        if ( hts_getline(in1, KS_SEP_LINE, &seq) <= 0 ) error("Truncated file? %s\n", args->fq1);
        if ( hts_getline(in1, KS_SEP_LINE, &str) <= 0 ) error("Truncated file? %s\n", args->fq1);
        if ( hts_getline(in1, KS_SEP_LINE, &qual) <= 0 ) error("Truncated file? %s\n", args->fq1);
        for (j=0; j<name.l; j++) if ( isspace(name.s[j]) ) { name.l = j; break; } 
        kputc('\n', &name);
        kputc('\n', &seq);
        kputc('\n', &str);
        kputc('\n', &qual);

        out_name.l = 0;
        kputs(name.s, &out_name);

        if ( ret >= 0 )  // keep the line, trim the sequence and add the BX:Z:ACGT... tag
        {
            char tmp[17];
            bxhash_key2seq(key, tmp);
            out_name.l--;
            kputs("\tBX:Z:", &out_name);
            kputs(tmp, &out_name);
            kputc('\n', &out_name);
            ssize_t ret = bgzf_write(out1->fp.bgzf, out_name.s, out_name.l);
            if ( ret!=out_name.l ) error("Failed to write %d bytes: %s\n", out_name.l,out1_fname);

            ret = bgzf_write(out1->fp.bgzf, seq.s + 16, seq.l - 16);
            if ( ret!=seq.l - 16 ) error("Failed to write %d bytes: %s\n", seq.l - 16,out1_fname);

            ret = bgzf_write(out1->fp.bgzf, str.s, str.l);
            if ( ret!=str.l ) error("Failed to write %d bytes: %s\n", str.l,out1_fname);

            ret = bgzf_write(out1->fp.bgzf, qual.s + 16, qual.l - 16);
            if ( ret!=qual.l - 16 ) error("Failed to write %d bytes: %s\n", qual.l - 16,out1_fname);

            nwr++;
            nwr_bx++;
            if ( ret==0 && args->reads.key[i] != key ) nwr_bx_mod++;
        }
        else if ( !args->discard )    // leave as is, no truncation
        {
            ret = bgzf_write(out1->fp.bgzf, out_name.s, out_name.l);
            if ( ret!=out_name.l ) error("Failed to write %d bytes: %s\n", out_name.l,out1_fname);
            ret = bgzf_write(out1->fp.bgzf, seq.s, seq.l);
            if ( ret!=seq.l ) error("Failed to write %d bytes: %s\n", seq.l,out1_fname);
            ret = bgzf_write(out1->fp.bgzf, str.s, str.l);
            if ( ret!=str.l ) error("Failed to write %d bytes: %s\n", str.l,out1_fname);
            ret = bgzf_write(out1->fp.bgzf, qual.s, qual.l);
            if ( ret!=qual.l ) error("Failed to write %d bytes: %s\n", qual.l,out1_fname);
            nwr++;
        }

        htsFile *in = in2 ? in2 : in1;
        if ( hts_getline(in, KS_SEP_LINE, &name2) <= 0 ) error("Truncated file? %s\n", in2 ? args->fq2 : args->fq1);
        if ( hts_getline(in, KS_SEP_LINE, &seq) <= 0 ) error("Truncated file? %s\n", in2 ? args->fq2 : args->fq1);
        if ( hts_getline(in, KS_SEP_LINE, &str) <= 0 ) error("Truncated file? %s\n", in2 ? args->fq2 : args->fq1);
        if ( hts_getline(in, KS_SEP_LINE, &qual) <= 0 ) error("Truncated file? %s\n", in2 ? args->fq2 : args->fq1);
        for (j=0; j<name2.l; j++) if ( isspace(name2.s[j]) ) { name2.l = j; break; }
        kputc('\n', &name2);
        kputc('\n', &seq);
        kputc('\n', &str);
        kputc('\n', &qual);
        // sanity check: the read names should match
        if ( strcmp(name.s,name2.s) ) error("Read pairs have different names? %s vs %s, line %d\n",name.s,name2.s,args->reads.n*nskip+1);
        if ( !args->discard || ret >= 0 )  // keep the line
        {
            ret = bgzf_write(out2->fp.bgzf, out_name.s, out_name.l);        // output the name of the first read, it has the BX tag attached
            if ( ret!=out_name.l ) error("Failed to write %d bytes: %s\n", out_name.l,out2_fname);
            ret = bgzf_write(out2->fp.bgzf, seq.s, seq.l);
            if ( ret!=seq.l ) error("Failed to write %d bytes: %s\n", seq.l,out2_fname);
            ret = bgzf_write(out2->fp.bgzf, str.s, str.l);
            if ( ret!=str.l ) error("Failed to write %d bytes: %s\n", str.l,out2_fname);
            ret = bgzf_write(out2->fp.bgzf, qual.s, qual.l);
            if ( ret!=qual.l ) error("Failed to write %d bytes: %s\n", qual.l,out2_fname);
        }
    }
    free(str.s);
    free(name.s);
    free(name2.s);
    free(out_name.s);
    free(seq.s);
    free(qual.s);
    if ( hts_close(in1)!=0 ) error("Close failed: %s\n", args->fq1);
    if ( in2 && hts_close(in2)!=0 ) error("Close failed: %s\n", args->fq2);
    if ( hts_close(out1)!=0 ) error("Close failed: %s\n", out1_fname);
    if ( hts_close(out2)!=0 ) error("Close failed: %s\n", out2_fname);
    free(out1_fname);
    free(out2_fname);

    ksprintf(&args->str,"%s.txt",args->prefix);
    args->fh_log = fopen(args->str.s,"w");
    if ( !args->fh_log ) error("Error opening %s for writing: %s\n",args->str.s,strerror(errno));

    fprintf(args->fh_log,"# nraw_reads              .. number of reads total\n");
    fprintf(args->fh_log,"# nraw_reads_known_bx     .. number of reads with a known barcode\n");
    fprintf(args->fh_log,"# nchanged_reads          .. number of reads with a modified barcode\n");
    fprintf(args->fh_log,"# ntrusted_reads          .. number of reads trusted because they come from a well populated barcode (%d reads)\n", args->trust_nreads);
    fprintf(args->fh_log,"# nfew_reads              .. number of reads that could be changed but were not because the barcode had fewer than %d reads\n", args->good_bx_nreads);
    fprintf(args->fh_log,"# nout_reads_known_bx     .. number of reads with a known barcode after error correction\n");
    fprintf(args->fh_log,"# nwr_reads_total         .. number of reads written to output total\n");
    fprintf(args->fh_log,"# nwr_reads_bx            .. number of reads written written with BX:Z:<barcode> added\n");
    fprintf(args->fh_log,"# nwr_reads_bx_changed    .. number of reads written written with BX:Z:<barcode> added and modified\n");
    fprintf(args->fh_log,"# trust_barcode_nreads    .. trust known barcodes with this many or more reads\n");
    fprintf(args->fh_log,"# good_barcode_nreads     .. barcodes changed only if the new barcode had at least this many reads\n");
    fprintf(args->fh_log,"# RAW_KNOWN               .. raw reads in known barcodes\n");
    fprintf(args->fh_log,"# RAW_UNKNOWN             .. raw reads in unknown barcode\n");
    fprintf(args->fh_log,"# OUT_KNOWN               .. reads in known barcodes after error correction\n");
    fprintf(args->fh_log,"# OUT_UNKNOWN             .. reads in unknown barcodes after error correction\n");
    fprintf(args->fh_log,"# MOD_NORI                .. the number of modified reads in the original barcodes\n");
    fprintf(args->fh_log,"# MOD_NNEW                .. the number of modified reads in the modified barcodes\n");
    fprintf(args->fh_log,"CMD\t%s", args->argv[0]);
    for (i=1; i<args->argc; i++) fprintf(args->fh_log," %s",args->argv[i]);
    fprintf(args->fh_log,"\n");

    uint64_t nraw_reads_known = 0, nout_reads_known = 0;
    uint32_t beg, end, n = dist_n(args->rpb_raw_kn);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->rpb_raw_kn, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "RAW_KNOWN\t%u\t%u\t%"PRIu64"\n", beg,end,cnt);   // nReads, not nBarcodes
        nraw_reads_known += cnt;
    }
    dist_destroy(args->rpb_raw_kn);

    n = dist_n(args->rpb_raw_un);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->rpb_raw_un, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "RAW_UNKNOWN\t%u\t%u\t%"PRIu64"\n", beg,end,cnt); // nReads, not nBarcodes
    }
    dist_destroy(args->rpb_raw_un);

    khint_t k;
    dist_t *dist = dist_init(4);
    for (k=0; k<kh_end(args->bx_out_kn); k++)
        if ( kh_exist(args->bx_out_kn,k) ) dist_insert(dist, kh_val(args->bx_out_kn,k));
    bxhash_destroy(args->bx_out_kn);
    n = dist_n(dist);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(dist, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "OUT_KNOWN\t%u\t%u\t%"PRIu64"\n", beg,end,(uint64_t)cnt*beg);
        nout_reads_known += cnt*beg;
    }
    dist_destroy(dist);

    dist = dist_init(4);
    for (k=0; k<kh_end(args->bx_out_un); k++)
        if ( kh_exist(args->bx_out_un,k) ) dist_insert(dist, kh_val(args->bx_out_un,k));
    bxhash_destroy(args->bx_out_un);
    n = dist_n(dist);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(dist, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "OUT_UNKNOWN\t%u\t%u\t%"PRIu64"\n", beg,end,(uint64_t)cnt*beg);
    }
    dist_destroy(dist);

    n = dist_n(args->rpb_mod_ori);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->rpb_mod_ori, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "MOD_NORI\t%u\t%u\t%"PRIu64"\n", beg,end,cnt); // nReads, not nBarcodes
    }
    dist_destroy(args->rpb_mod_ori);

    n = dist_n(args->rpb_mod_new);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->rpb_mod_new, i, &beg, &end);
        if ( !cnt ) continue;
        fprintf(args->fh_log, "MOD_NNEW\t%u\t%u\t%"PRIu64"\n", beg,end,cnt); // nReads, not nBarcodes
    }
    dist_destroy(args->rpb_mod_new);

    fprintf(args->fh_log,"SN\tnraw_reads\t%"PRIu64"\n", args->reads.n);
    fprintf(args->fh_log,"SN\tnraw_reads_known_bx\t%"PRIu64"\t%.2f%%\n", nraw_reads_known, args->reads.n ? nraw_reads_known*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnchanged_reads\t%"PRIu64"\t%.2f%%\n", nchanged, args->reads.n ? nchanged*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tntrusted_reads\t%"PRIu64"\t%.2f%%\n", ntrusted, args->reads.n ? ntrusted*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnfew_reads\t%"PRIu64"\t%.2f%%\n", args->nfew_reads, args->reads.n ? args->nfew_reads*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnout_reads_known_bx\t%"PRIu64"\t%.2f%%\n", nout_reads_known, args->reads.n ? nout_reads_known*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnwr_reads_total\t%"PRIu64"\t%.2f%%\n", nwr, args->reads.n ? nwr*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnwr_reads_bx\t%"PRIu64"\t%.2f%%\n", nwr_bx, args->reads.n ? nwr_bx*100./args->reads.n : 0);
    fprintf(args->fh_log,"SN\tnwr_reads_bx_changed\t%"PRIu64"\t%.2f%%\n", nwr_bx_mod, args->reads.n ? nwr_bx_mod*100./args->reads.n : 0);
    fprintf(args->fh_log,"LM\ttrust_barcode_nreads\t%d\n", args->trust_nreads);
    fprintf(args->fh_log,"LM\tgood_barcode_nreads\t%d\n", args->good_bx_nreads);
}

static void usage(void)
{
    fprintf(stderr,"\n");
    fprintf(stderr,"About:  Extract barcode sequences and correct sequencing errors\n");
    fprintf(stderr,"Usage:  bxcheck trim [OPTIONS] [in.12.fq.gz|in.1.fq.gz in.2.fq.gz]\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"    -D, --dont-discard              Do not discard unassignable reads\n");
    fprintf(stderr,"    -g, --good-barcodes INT         Switch only to well-populated barcodes with at least INT reads [5]\n");
    fprintf(stderr,"    -l, --barcodes-list FILE        List of 10x barcodes\n");
    fprintf(stderr,"    -o, --output <prefix>           Write output to prefix.1.fq.gz and prefix.2.fq.gz\n");
    fprintf(stderr,"    -t, --trust-barcodes INT        Trust barcodes with this many or more reads [40]\n");
    fprintf(stderr,"\n");
    exit(1);
}

int main_trim(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc = argc;
    args->argv = argv;
    args->trust_nreads = 40;    // For speed up: 99.78% of modified reads are from barcodes with <=40 reads, save 88% of comparisons
    args->good_bx_nreads = 5;
    args->discard = 1;
    base_qual_init_qbin2prob();

    static const struct option loptions[] =
    {
        {"dont-discard", required_argument, NULL, 'D'},
        {"trust-barcodes", required_argument, NULL, 't'},
        {"good-barcodes", required_argument, NULL, 'g'},
        {"output", required_argument, NULL, 'o'},
        {"barcodes-list", required_argument, NULL, 'l'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ( (opt=getopt_long(argc,argv,"?hl:o:t:Dg:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case 'D': args->discard = 0; break;
            case 't': args->trust_nreads = atoi(optarg); break;
            case 'g': args->good_bx_nreads = atoi(optarg); break;
            case 'o': args->prefix = optarg; break;
            case 'l': args->bx_known = bxhash_init(optarg); break;
            case '?':
            case 'h': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fq1 = "-";  // reading from stdin
        else usage();
    }
    else
    {
        args->fq1 = argv[optind];
        if ( optind+1<argc ) args->fq2 = argv[optind+1];
    }

    if ( !args->bx_known ) error("Missing the -l, --barcodes-list option\n");
    if ( !args->prefix ) error("Missing the -o, --output option\n");

    read_fq_bx(args);
    analyze_fq_bx(args);

    destroy(args);
    return 0;
}

