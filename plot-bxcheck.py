#!/usr/bin/env python
#
#   Author: petr.danecek@sanger
#   About:  Script for plotting the output of bxcheck
#   Usage:  
#       bxcheck file.bam > file.txt
#       plot-bxcheck file.txt -d dir
#

import sys

def usage():
    print 'Usage: plot-bxcheck.py [OPTIONS] bx-output.txt'
    print 'Options:'
    print '   -d, --dir STR     output directory'
    print '   -n, --name STR    main title [bxcheck stats]'
    print '   -t, --type STR    image type [png]'
    sys.exit(1)


fname = None
dir   = None
itype = 'png'
main_title = 'bxcheck stats'
replot = True

if len(sys.argv) < 2: usage()
args = sys.argv[1:]
while len(args):
    if args[0]=='-d' or args[0]=='--dir': 
        args = args[1:]
        dir  = args[0]
    elif args[0]=='-t' or args[0]=='--type': 
        args  = args[1:]
        itype = args[0]
    elif args[0]=='-n' or args[0]=='--name': 
        args  = args[1:]
        main_title = args[0]
    elif args[0]=='--dont-replot':  # for debugging, do not call matplotlib if data file older than existing image
        replot = False
    else:
        fname = args[0]
    args = args[1:]

if dir==None: usage()

color = [
    '#d43f3a',   # coral red
    '#06abb4',   # royal blue

#    '#a7c65c',   # aloe
#    '#ee4789',   # pink
    
#    '#a7c65c',   # aloe
#    '#2e5251',   # pine
]

graphs = []

dists = \
{
    'FRAG_NREADS': 
    {
        'ylabel': 'Cumulative frequency',
        'xlabel': 'Reads per fragment',
        'img':    'frag_nreads',
        'hist':   1,
        'xlim':   0.99,
        'plot_density': 1,
        'cfrac':  1,
    },
    'FRAG_SIZE': 
    {
        'ylabel': 'Cumulative frequency',
        'xlabel': 'Fragment length',
        'img':    'frag_size',
        'xlim':   0.90,
        'hist':   1,
        'xsci':   1,
        'plot_density': 1,
        'cfrac':  1,
    },
}

dists2 = \
{
    'DIST_ALL_SCLIP1': { 'img':'dist_sclip', 'dat':'dat1' },
    'DIST_ALL_SCLIP2': { 'img':'dist_sclip', 'dat':'dat2' },
    'dist_sclip':
    {
        'label1': 'First reads',
        'label2': 'Second reads',
        'xlabel': 'Number of bases after soft-clipping',
        'ylabel': 'Number of reads',
        'ylog':   1,
    },
    'DIST_COV_ALL_READS':  { 'img':'cov_reads', 'dat':'dat1' },
    'DIST_COV_GOOD_READS': { 'img':'cov_reads', 'dat':'dat2' },
    'cov_reads':
    {
        'ylabel': 'Number of sites',
        'xlabel': 'Coverage',
        'label1': 'All mapped reads',
        'label2': 'Good reads',
        'hist':   1,
        'xlim':   0.99,
        'ysci':   1,
        'alpha': 0.8,
    },
    'BX_NALL_READS':  { 'img':'bx_reads', 'dat':'dat1' },
    'BX_NGOOD_READS': { 'img':'bx_reads', 'dat':'dat2' },
    'bx_reads':
    {
        'ylabel': 'Cumulative frequency',
        'xlabel': 'Reads per barcode',
        'label1': 'All reads',
        'label2': 'Good reads',
        'xlim':   0.92,
        'cfrac':  1,
        'alpha': 0.8,
        'plot_density': 1,
    },
}


import matplotlib as mpl
from matplotlib.image import BboxImage
from matplotlib.transforms import Bbox, TransformedBbox
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv, os, numpy, base64
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

try:
    os.makedirs(dir)
except:
    pass


def get_xlim(xdat,ydat,frac):
    xlim = 0
    sum = 0
    tot = 0
    for cnt in ydat: tot += cnt
    for i in range(len(xdat)):
        sum += ydat[i]
        xlim = xdat[i]
        if sum/tot > frac: break
    return xlim

def plot_dist(dist):
    xdat = dist['xdat']
    ydat = dist['ydat']
    
    wh = (7,3)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')
 
    if 'cfrac' in dist:
        xlim = 1
        if 'xlim' in dist: xlim = dist['xlim']
        frac = 0
        tot  = sum(ydat)
        i = 0
        for i in range(len(xdat)):
            frac += ydat[i]/tot
            ydat[i] = frac
            if frac >= xlim: break
        xdat = xdat[0:i]
        ydat = ydat[0:i]
        ax1.plot(xdat,ydat,'o-',color=color[0],mec=color[0],alpha=0.8)
    else:
        width = xdat[1] - xdat[0]
        if len(xdat) > 100:
            ax1.fill_between([x-width*0.5 for x in xdat],[0 for y in ydat],ydat,color=color[0])
        else:
            ax1.bar([x-width*0.5 for x in xdat],ydat,width=width,color=color[0],edgecolor=color[0])
        if 'xlim' in dist:
            xlim = get_xlim(xdat,ydat,dist['xlim'])
            ax1.set_xlim(-width,xlim)
        else:
            ax1.set_xlim(-width)
    
    
    if 'xsci' in dist and dist['xsci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='x')
    
    if 'ysci' in dist and dist['ysci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    
    if 'ylog' in dist:
        ax1.set_yscale('symlog')

    if 'ylabel' in dist: ax1.set_ylabel(dist['ylabel'])
    if 'xlabel' in dist: ax1.set_xlabel(dist['xlabel'])

    plt.subplots_adjust(bottom=0.15)
    plt.savefig(dir+'/'+dist['img']+'.'+itype)
    plt.close()

def dat_to_cfrac(xdat,ydat,xlim):
    frac = 0
    tot  = sum(ydat)
    i = 0
    for i in range(len(xdat)):
        if xdat[i] > xlim: break
        frac += ydat[i]/tot
        ydat[i] = frac
    xdat = xdat[0:i]
    ydat = ydat[0:i]
    return (xdat,ydat)

def plot_dist2(dist):
    wh = (7,3)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')

    xdat1 = dist['dat1']['xdat']
    ydat1 = dist['dat1']['ydat']
    xdat2 = dist['dat2']['xdat']
    ydat2 = dist['dat2']['ydat']

    line1 = 'o-'
    line2 = 'o-'
    
    col2 = color[0]
    col1 = color[1]
    
    alpha = 0.8
    if 'alpha' in dist: alpha = dist['alpha']

    plots  = None
    labels = [dist['label1'],dist['label2']]
    loc = 'best'

    if 'hist' in dist and dist['hist']==1:
        width = xdat1[1] - xdat1[0]
        if width > xdat2[1] - xdat2[0]: width = xdat2[1] - xdat2[0]
        if len(xdat1) > 100 or len(xdat2) > 100:
            ax1.fill_between([x-width*0.5 for x in xdat1],[0 for y in ydat1],ydat1,color=col1,alpha=alpha)
            ax1.fill_between([x-width*0.5 for x in xdat2],[0 for y in ydat2],ydat2,color=col2,alpha=alpha)
            plots = [plt.Rectangle((0,0),1,1,color=col1,alpha=alpha),plt.Rectangle((0,0),1,1,color=col2,alpha=alpha)]
        else:
            plots  = ax1.bar([x-width*0.5 for x in xdat1],ydat1,width=width,color=col1,edgecolor=col1,alpha=alpha)
            plots += ax1.bar([x-width*0.5 for x in xdat2],ydat2,width=width,color=col2,edgecolor=col2,alpha=alpha)
        if 'xlim' in dist:
            xlim1 = get_xlim(xdat1,ydat1,dist['xlim'])
            xlim2 = get_xlim(xdat2,ydat2,dist['xlim'])
            if xlim1 < xlim2: xlim1 = xlim2
            ax1.set_xlim(-width,xlim1)
        else:
            ax1.set_xlim(-width)
    elif 'cfrac' in dist and dist['cfrac']==1:
        loc = 'lower right'
        xlim = 1
        if 'xlim' in dist: 
            xlim  = get_xlim(xdat1,ydat1,dist['xlim'])
            xlim2 = get_xlim(xdat2,ydat2,dist['xlim'])
            if xlim < xlim2: xlim = xlim2
        (xdat1,ydat1) = dat_to_cfrac(xdat1,ydat1,xlim)
        (xdat2,ydat2) = dat_to_cfrac(xdat2,ydat2,xlim)
        plots  = ax1.plot(xdat1,ydat1,line1,color=col1,mec=col1,alpha=alpha)
        plots += ax1.plot(xdat2,ydat2,line2,color=col2,mec=col2,alpha=alpha)
    else:
        plots  = ax1.plot(dist['dat1']['xdat'],dist['dat1']['ydat'],line1,color=col1,mec=col1,alpha=alpha)
        plots += ax1.plot(dist['dat2']['xdat'],dist['dat2']['ydat'],line2,color=col2,mec=col2,alpha=alpha)
    
    plt.legend(plots,labels,numpoints=1,markerscale=1,loc=loc,prop={'size':10},frameon=False)
    if 'ylog' in dist and dist['ylog']==1: ax1.set_yscale('symlog')
    ax1.set_ylabel(dist['ylabel'])
    ax1.set_xlabel(dist['xlabel'])
    if 'title' in dist:
        ax1.set_title(dist['title'])
    if 'ysci' in dist and dist['ysci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    plt.subplots_adjust(bottom=0.15)
    if itype!='png': plt.savefig(dir+'/'+dist['img']+'.'+itype)
    plt.savefig(dir+'/'+dist['img']+'.png')
    plt.close()


def parse_file(fname):
    if fname==None: f = sys.stdin
    else: f = open(fname, 'r')
    dat = {}
    dat['SN'] = {}
    dat['LM'] = {}
    reader = csv.reader(f, 'tab')
    for row in reader:
        for key in dists:
            if row[0] == key:
                if key not in dat: dat[key] = {}
                if 'xdat' not in dat[key]:
                    dat[key]['xdat'] = []
                    dat[key]['ydat'] = []
                    dat[key]['rdat'] = []
                    dat[key].update(dists[key])
                dat[key]['xdat'].append(float(row[1]))
                dat[key]['rdat'].append(float(row[3]))
                if 'plot_density' in dists[key]:
                    dat[key]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
                else:
                    dat[key]['ydat'].append(float(row[3]))
        if row[0]=="SN":
            dat['SN'][row[1]] = row[2]
            if len(row)==4: dat['LM'][row[1]] = row[3]
        if row[0]=="LM":
            dat['LM'][row[1]] = row[2]
        for key in dists2:
            if row[0] == key:
                img  = dists2[key]['img']
                key2 = dists2[key]['dat']
                if img not in dat: dat[img] = {}
                if key2 not in dat[img]: 
                    dat[img][key2] = {}
                    dat[img][key2]['ydat'] = []
                    dat[img][key2]['xdat'] = []
                    dat[img].update(dists2[img])
                dat[img][key2]['xdat'].append(float(row[1]))
                if 'plot_density' in dat[img]:
                    dat[img][key2]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
                else:
                    dat[img][key2]['ydat'].append(float(row[2])) 
                dat[img]['img'] = img
    return dat

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

def percent(part, total):
    return "%.1f%%" % (float(part)*100./float(total));

def embed_image(name):
    fh = open(dir+'/'+name+'.'+itype, "rb")
    return "<img style='width:100%;margin-top:1em;margin-bottom:1em;' src='data:image/png;base64," + base64.b64encode(fh.read()) + "'>"

help_id = 0
def help_text(text):
    global help_id
    help_id = help_id + 1
    return \
        "<div style='float:right;padding:0.5em;'><div onclick='toggle_help(event,"+str(help_id)+")' class='help'>?</div></div>" + \
        "<div onclick='toggle_help(event,"+str(help_id)+")' id='help"+str(help_id)+"' class='help_text'>"+text+"</div>"

def write_html(fname, dists):
    excluded = ''
    excluded_help = ''
    if 'n_excluded_soft_clips' in dat['SN']:
        excluded = """<tr><td class='lalni'>.. soft-clips&gt;"""+dat['LM']['n_excluded_soft_clips'] + """<td>""" \
                            + bignum(dat['SN']['n_excluded_soft_clips']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_soft_clips'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help = """
            <dt>Soft-clips</dt>
            <dd>Number of reads excluded because of too many soft-clipped bases.</dd>"""
    if 'n_excluded_bx_has_n' in dat['SN'] and dat['SN']['n_excluded_bx_has_n']!=0:
        excluded = """<tr><td class='lalni'>.. barcode has N<td>""" \
                            + bignum(dat['SN']['n_excluded_bx_has_n']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_bx_has_n'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help = """
            <dt>barcode has N</dt>
            <dd>Unknown base in the barcode sequence.</dd>"""

    filters = ''
    filters_help = ''
    if 'n_flag_DUP' in dat['LM']:
        filters += """<tr><td class='lalnii'>.. duplicates<td>""" \
                    + bignum(dat['SN']['n_flag_DUP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_DUP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Duplicates</dt>
            <dd>Number of reads marked as duplicate.</dd>"""
    if 'n_flag_UNMAP' in dat['LM']:
        filters += """<tr><td class='lalnii'>.. unmapped<td>""" \
                    + bignum(dat['SN']['n_flag_UNMAP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_UNMAP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Unmapped</dt>
            <dd>Number of unmapped reads.</dd>"""
    if 'n_flag_MUNMAP' in dat['LM']:
        filters += """<tr><td class='lalnii'>.. mate unmapped<td>""" \
                    + bignum(dat['SN']['n_flag_MUNMAP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_MUNMAP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Mate unmapped</dt>
            <dd>Number of reads with unmapped mates.</dd>"""
    if 'n_flag_SUPPLEMENTARY' in dat['LM']:
        filters += """<tr><td class='lalnii'>.. supplementary<td>""" \
                    + bignum(dat['SN']['n_flag_SUPPLEMENTARY']) + """<td>""" \
                    + percent(dat['SN']['n_flag_SUPPLEMENTARY'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Supplementary</dt>
            <dd>Number of supplementary reads.</dd>"""
    if 'n_flag_SECONDARY' in dat['LM'] and dat['SN']['n_flag_SECONDARY']!='0':
        filters += """<tr><td class='lalnii'>.. secondary<td>""" \
                    + bignum(dat['SN']['n_flag_SECONDARY']) + """<td>""" \
                    + percent(dat['SN']['n_flag_SECONDARY'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Secondary</dt>
            <dd>Number of secondary reads.</dd>"""
    if 'n_flag_QCFAIL' in dat['LM'] and dat['SN']['n_flag_QCFAIL']!='0':
        filters += """<tr><td class='lalnii'>.. QC failed<td>""" \
                    + bignum(dat['SN']['n_flag_QCFAIL']) + """<td>""" \
                    + percent(dat['SN']['n_flag_QCFAIL'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>QC failed</dt>
            <dd>Number of reads not passing quality controls.</dd>"""

    fh = open(fname,"w")
    fh.write("""<!DOCTYPE html><html>
        <style>
            .column {
                clear: left;
                width: 25%;
                display: inline-block;
                vertical-align: top;
            }
            .box {
                padding: 0.5em;
                margin: 0.5em;
                margin-bottom: 1em;
                border: solid 1px #ddd;
                border-radius: 0.2em;
            }
            h2 {
            }
            dt {
                font-weight: bold;
            }
            dd {
                margin-bottom: 1em;
                margin-left: 1em;
            }
            .help {
                display: inline-block;
                margin-left:0.5em;
                width: 1.5em;
                height: 1.5em;
                line-height: 1.5em;
                font-size: 0.5em;
                font-weight: bold;
                text-align: center;
                color:white;
                background-color: #ddd;
                border-radius: 50%;
                cursor: pointer;
            }
            .help_text {
                font-size: small;
                font-weight: normal;
                text-align: left;
                width: 90%;
                display: none;
                background-color: white;
                border-radius: 0.2em;
                border: solid 1px #ddd;
                cursor: pointer;
                padding: 1em;
                box-shadow: 3px 3px 3px #ddd;
            }
            .laln { text-align: left; }
            .raln { text-align: right; }
            .lalni { text-align: left; padding-left: 1.5em; }
            .lalnii { text-align: left; padding-left: 3em; }
            table.numbers {
                width: 90%;
                margin: auto;
                text-align: right;
            }
            div.sep {
                margin:2em;
                margin-left: auto;
                margin-right: auto;
                width: 80%; 
                height: 1px; 
                background: #ddd;
                overflow: hidden;
            }
        </style>
        <script type="text/javascript">
            function toggle_help(e, id)
            {
                div = document.getElementById('help'+id)
                if ( div.style.display=='block' ) div.style.display='none';
                else div.style.display='block';
            }
        </script>
        <body> 
            <div style="width:100%; text-align:center;">
            <span class="column">
                <div class="box">
                    """+ help_text("""
                        <dl>
                            <dt>Reads total</dt>
                                <dd>Total number of reads, each read pair counts as two.</dd>
                            <dt>Barcoded reads</dt>
                                <dd>Reads with the BX tag</dd>
                            <dt>Reads excluded</dt>
                                <dd>Reads excluded because of flag, low mapping quality, etc.</dd>
                            <dt>MQ&lt;"""+dat['LM']['n_excluded_mq'] + """</dt>
                                <dd>Number of reads excluded because of low mapping quality.</dd>
                            """ + excluded_help + """
                            <dt>Filtered</dt>
                                <dd>Number of reads excluded because of BAM flags.</dd>
                            """ + filters_help + """
                            <dt>Coverage</dt>
                                <dd>Good reads were filtered by flag (""" + dat['LM']['flags'] + """) and mapping quality (&ge;""" + dat['LM']['minMQ'] + """).
                                    The x axis range was set to show at least """ + str(100*dat['cov_reads']['xlim']) + """% of the data.
                                </dd>
                            <dt>Number of bases after soft-clipping</dt>
                                <dd>Histogram of read-lengths after trimming soft-clipped bases. Includes all reads. </dd>
                        </dl>
                    """) +"""
                    <h2>Sequencing</h2>
                    <table class='numbers'>
                        <tr><td class='laln'>Reads total <td>""" 
                            + bignum(dat['SN']['n_all_reads']) + """<td></tr>
                        <tr><td class='laln'>Barcoded reads <td>"""
                            + bignum(dat['SN']['n_bx_reads']) + """<td>"""
                            + percent(dat['SN']['n_bx_reads'],dat['SN']['n_all_reads']) +"""</tr>
                        <tr><td class='laln'>Reads excluded <td>"""
                            + bignum(dat['SN']['n_excluded']) + """<td>""" 
                            + percent(dat['SN']['n_excluded'],dat['SN']['n_all_reads']) + """</tr>
                        <tr><td class='lalni'>.. MQ&lt;"""+dat['LM']['n_excluded_mq'] + """<td>"""
                            + bignum(dat['SN']['n_excluded_mq']) + """<td>""" 
                            + percent(dat['SN']['n_excluded_mq'],dat['SN']['n_all_reads']) + """</tr>
                        """ + excluded + """
                        <tr><td class='lalni'>.. filtered<td>"""
                            + bignum(dat['SN']['n_excluded_flag']) + """<td>""" 
                            + percent(dat['SN']['n_excluded_flag'],dat['SN']['n_all_reads']) + """</tr>
                        """ + filters + """
                    </table>
                    <div class='sep'></div> """
                    + embed_image(dat['cov_reads']['img']) + """ <div class='sep'></div> """
                    + embed_image(dat['dist_sclip']['img']) + """</div>
            </span>
            <span class="column">
                <div class="box">
                    """+ help_text("""
                        <dl>
                            <dt>GEM</dt>
                                <dd>Gel Bead in Emulsion. In the sequencing data GEM corresponds to a set of reads with the same BX tag.</dd>
                            <dt>Excluded</dt>
                                <dd>Number of excluded barcodes</dd>
                            <dt>Reads per barcode</dt>
                                <dd>Cumulative frequency of barcodes with this many or fewer
                                    reads per barcode. The x axis range was set to include """ + str(100*dat['bx_reads']['xlim']) + """% of the data.
                                </dd>
                            <dt>Reads per fragment</dt>
                                <dd>Cumulative frequency of fragments with this many or fewer reads.
                                The x axis range was set to include """ + str(100*dat['FRAG_NREADS']['xlim']) + """% of the data.
                                </dd>
                        </dl>
                    """) +"""
                    <h2>GEMs</h2>
                    <table class='numbers'>
                        <tr><td class='laln'>GEMs total <td>""" 
                            + bignum(dat['SN']['n_all_barcodes']) + """<td><td></tr>
                        <tr><td class='laln'>Excluded<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='lalni'>.. fewer than """ + dat['LM']['n_barcodes_excluded_min_reads'] + """ reads in the barcode<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded_min_reads']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded_min_reads'],dat['SN']['n_all_barcodes']) +"""</tr>
                    </table>
                    <div class='sep'></div> """
                    + embed_image(dat['bx_reads']['img']) + """ <div class='sep'></div> """
                    + embed_image(dat['FRAG_NREADS']['img']) + """</div>
                <div class="box">
                    """+ help_text("""
                        <dl>
                            <dt>Molecule length</dt>
                                <dd>Cumulative frequency of estimated DNA fragment size. The fragment
                                size is estimated simply as the distance between the first and last
                                reads within the barcode that to the same chromosome.
                                The x axis range was set to include """ + str(100*dat['FRAG_SIZE']['xlim']) + """% of the data.
                                </dd>
                        </dl>
                    """)
                    + """<h2>Input DNA</h2>""" 
                    + embed_image(dat['FRAG_SIZE']['img']) + """</div>
            </span>
            </div>
        </body></html>
        """)
    fh.close()

#--------------------

plots = []
dat  = parse_file(fname)

for key in dists2:
    if key not in dat: continue
    if 'img' in dists2[key]: continue
    if key not in graphs: graphs.append(key)
    if not replot:
        img = dir+'/'+dat[key]['img']+'.'+itype
        if os.path.exists(img) and os.path.getctime(img) > os.path.getctime(fname): continue
    plot_dist2(dat[key])

for key in dists:
    if key not in dat: continue
    if key not in graphs: graphs.append(key)
    if not replot:
        img = dir+'/'+dat[key]['img']+'.'+itype
        if os.path.exists(img) and os.path.getctime(img) > os.path.getctime(fname): continue
    plot_dist(dat[key])

#   for key in chr_maps:
#       if key not in dat: continue
#       if key not in graphs: graphs.append(key)
#       plot_dist(dat[key])
#   
#   for key in graphs:
#       if key not in dat: continue
#       plots.append(dat[key])

write_html(dir+'/'+"plots.html", plots)



