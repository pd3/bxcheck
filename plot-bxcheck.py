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
    else:
        fname = args[0]
    args = args[1:]

if dir==None: usage()

graphs = ['dist_sclip','DIST_BX_NREADS','DIST_NFRAGS','DIST_FRAG_SIZE','DIST_FRAG_NREADS','DIST_READ_SPACING','COV_READS','COV_FRAGS']

dists = \
{
    'DIST_BX_NREADS': 
    {
        'title':  'Number of reads per barcode',
        'xlabel': 'Reads per barcode',
        'img':    'dist_bx_nreads',
        'xlog':   1,
    },
    'DIST_NFRAGS': 
    {
        'xlabel': 'Fragments per barcode',
        'title':  'Number of fragments per barcode',
        'img':    'dist_nfrags',
        'xlog':   1,
    },
    'DIST_FRAG_NREADS': 
    {
        'xlabel': 'Reads per fragment',
        'title':  'Number of reads per fragment',
        'img':    'dist_frag_nreads',
        'xlog':   1,
    },
    # 'DIST_READ_SPACING': 
    # {
    #     'xlabel': 'Read spacing within fragment',
    #     'title':  'Read spacing within a fragment',
    #     'img':    'dist_min_read_spacing',
    #     'xlog':   1,
    # },
    'DIST_FRAG_SIZE': 
    {
        'xlabel': 'Fragment size',
        'title':  'Fragment size distribution',
        'img':    'dist_frag_size',
        'xlog':   1,
    },
}

dists4 = \
{
    'DIST_SCLIP1': { 'img':'dist_sclip', 'dat':'dat3' },
    'DIST_SCLIP2': { 'img':'dist_sclip', 'dat':'dat4' },
    'DIST_ALL_SCLIP1': { 'img':'dist_sclip', 'dat':'dat1' },
    'DIST_ALL_SCLIP2': { 'img':'dist_sclip', 'dat':'dat2' },
    'dist_sclip':
    {
        'title': 'Read length distribution after soft-clipping',
        'img': 'dist_sclip',
        'label1': 'First reads',
        'label2': 'Second reads',
        'label3': 'First reads (filtered)',
        'label4': 'Second reads (filtered)',
        'xlabel': 'Number of bases after soft-clipping',
        'ylabel': 'Fraction of reads',
    },
}

chr_maps = \
{
    'COV_READS': 
    {
        'xlabel': 'Reads per %.0f kb',
        'ylabel': 'Frequency',
        'title':  'Genome coverage (reads)',
        'img': 'dist_cov_reads',
        'chr': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'],
        'xlog':   1,
    },
    'COV_FRAGS': 
    {
        'xlabel': 'Fragments per %.0f kb',
        'ylabel': 'Frequency',
        'title':  'Genome coverage (fragments)',
        'img': 'dist_cov_frags',
        'chr': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'],
        'xlog':   1,
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

def plot_dist(dist):
    xdat = dist['xdat']
    ydat = dist['ydat']
    cdat = []
    sum  = 0
    norm = 0
    for y in dist['rdat']:
        norm += float(y)
    for y in dist['rdat']:
        sum += float(y)
        cdat.append(sum/norm)
    
    wh = (7,5)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax2 = ax1.twinx()
    
    label1 = 'Cumulative fraction'
    label2 = 'Density'
    
    line1 = '.-'
    line2 = '.-'
    
    col1 = 'black'    # +col1 #4CAE4C (green), #EEA236 (orange), #D43F3A (red), #2E6DA4 (blue)
    col2 = '#D43F3A'
    
    p1 = ax1.plot(xdat,cdat,line2,label=label1,color=col1)
    p2 = ax2.plot(xdat,ydat,line1,label=label2,color=col2)
    plots = p1 + p2
    
    ax2.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    if 'xsci' in dist and dist['xsci']:
        ax2.ticklabel_format(style='sci', scilimits=(-2,2), axis='x')
    
    ax1.set_ylabel(label1, color=col1)
    ax2.set_ylabel(label2, color=col2)
    
    for tl in ax1.get_yticklabels(): tl.set_color(col1)
    for tl in ax2.get_yticklabels(): tl.set_color(col2)
    
    xlabel = dist['xlabel']
    ax1.set_xlabel(xlabel)
    
    ax2.set_yscale('log')

    if 'xlog' in dist and dist['xlog']:
        ax1.set_xscale('symlog')
        ax2.set_xscale('symlog')
    
    title = dist['title']
    ax1.set_title(title)
    
    labels = [l.get_label() for l in plots]
    plt.legend(plots,labels,numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
    
    # plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95
    
    plt.savefig(dir+'/'+dist['img']+'.'+itype)
    plt.close()


def plot_dist4(dist):
    wh = (7,5)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    
    label1 = dist['label1']
    label2 = dist['label2']
    label3 = dist['label3']
    label4 = dist['label4']
    
    line1 = '.-'
    line2 = '.-'
    line3 = '.-'
    line4 = '.-'
    
    col1 = 'black'    # +col1 #4CAE4C (green), #EEA236 (orange), #D43F3A (red), #2E6DA4 (blue)
    col2 = '#D43F3A'
    col3 = '#4CAE4C'
    col4 = '#EEA236'
    
    ax1.plot(dist['dat1']['xdat'],dist['dat1']['ydat'],line1,label=label1,color=col1)
    ax1.plot(dist['dat2']['xdat'],dist['dat2']['ydat'],line2,label=label2,color=col2)
    ax1.plot(dist['dat3']['xdat'],dist['dat3']['ydat'],line3,label=label3,color=col3)
    ax1.plot(dist['dat4']['xdat'],dist['dat4']['ydat'],line4,label=label4,color=col4)
    
    ax1.set_ylabel(dist['ylabel'])
    ax1.set_xlabel(dist['xlabel'])
    ax1.set_title(dist['title'])
    plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)
    
    # plt.subplots_adjust(**adjust)        # for example: bottom=0.2,left=0.1,right=0.95
    
    if itype!='png': plt.savefig(dir+'/'+dist['img']+'.'+itype)
    plt.savefig(dir+'/'+dist['img']+'.png')
    plt.close()


def parse_file(fname):
    dat = {}
    dat['SN'] = {}
    dat['LM'] = {}
    with open(fname, 'rb') as f:
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
                    dat[key]['xdat'].append(row[1])
                    dat[key]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
                    dat[key]['rdat'].append(row[3])
            if row[0]=="SN":
                dat['SN'][row[1]] = row[2]
                if len(row)==4: dat['LM'][row[1]] = row[3]
            for key in dists4:
                if row[0] == key:
                    img  = dists4[key]['img']
                    key2 = dists4[key]['dat']
                    if img not in dat: dat[img] = {}
                    if key2 not in dat[img]: 
                        dat[img][key2] = {}
                        dat[img][key2]['ydat'] = []
                        dat[img][key2]['xdat'] = []
                        dat[img].update(dists4[img])
                    dat[img][key2]['xdat'].append(row[1])
                    dat[img][key2]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
            for key in chr_maps:
                if row[0] != key: continue
                chr = row[1]
                beg = float(row[2])
                end = float(row[3])
                cnt = int(row[4])
                if chr not in chr_maps[key]['chr']: continue
                if key not in dat: 
                    dat[key] = {}
                    dat[key]['xdat'] = []
                    dat[key]['ydat'] = []
                    dat[key]['rdat'] = []
                if 'win' not in dat[key]: 
                    dat[key]['win']  = float(end) - float(beg)
                    dat[key]['prev'] = {'chr':chr,'end':0}
                if chr != dat[key]['prev']['chr']:
                    dat[key]['prev']['chr'] = chr
                    dat[key]['prev']['end'] = 0
                if key not in dat: 
                    dat[key] = {}
                if 'dist' not in dat[key]:
                    dat[key]['dist'] = { 0:0 }
                if 'map' not in dat[key]:
                    dat[key]['map'] = {}
                    if chr in dat[key]['map']: dat[key]['map'][chr] = []
                while dat[key]['prev']['chr']==chr and dat[key]['prev']['end'] < beg:
                    dat[key]['dist'][0] += 1
                    dat[key]['prev']['end'] += dat[key]['win']

                if cnt not in dat[key]['dist']:
                    dat[key]['dist'][cnt] = 1
                else:
                    dat[key]['dist'][cnt] += 1

    # Plot soft clips as cumulative distribution of read length after clipping
    for key in ['dat1','dat2','dat3','dat4']:
        norm = 0.0
        for y in dat['dist_sclip'][key]['ydat']: norm += y
        sum = 0
        for i in range(len(dat['dist_sclip'][key]['xdat'])):
            sum += dat['dist_sclip'][key]['ydat'][i]
            dat['dist_sclip'][key]['ydat'][i] = sum / norm

    # total number of fragments
    dat['SN']['n_fragments'] = 0
    for i in range(len(dat['DIST_NFRAGS']['rdat'])):
        dat['SN']['n_fragments'] += int(dat['DIST_NFRAGS']['rdat'][i]) * int(dat['DIST_NFRAGS']['xdat'][i])

    # depth distribution
    for key in chr_maps:
        if key not in dat: continue
        dat[key].update(chr_maps[key])
        dat[key]['xlabel'] = dat[key]['xlabel'] % (dat[key]['win']*1e-3)
        dat[key]['xdat'] = []
        dat[key]['ydat'] = []
        dat[key]['rdat'] = []
        out = sorted(dat[key]['dist'].keys())
        for cnt in out:
            dat[key]['xdat'].append(cnt)
            dat[key]['ydat'].append(dat[key]['dist'][cnt])
            dat[key]['rdat'].append(dat[key]['dist'][cnt])

    return dat

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

def percent(part, total):
    return "%.1f%%" % (float(part)*100./float(total));

def write_html(fname, dists):
    fh = open(fname,"w")
    fh.write("""<!DOCTYPE html><html>
        <style>
            li {
               margin: 0 0 0.5em 0;
            }
            table.summary {
                position: absolute;
                top:2em;
                left:5em;
            }
            table.summary td {
                padding: 0.2em 1em 0.2em 1em;
                text-align:left;
            }
            table.summary td.indent {
                padding-left:4em;
            }
        </style>
        <body> 
        <div style="width:100%; text-align:center;">
        <table style="display:inline-block;">
        <tr><td colspan="2" style="background-color:#DFF0D8; padding:0.2em;">
                <h3>""" + main_title + """</h3>
            </td>
        <tr>
            <td style="background-color:#FCF8E3;vertical-align:top;padding:1em;text-align:left;position:relative;">

            <p><b>Stats</b>
            <ul>
            <li><span style='cursor:pointer' onclick='show_div(0)'>Summary numbers</span>
            </ul>

            <p><b>Graphs</b>
            <ul>
        """)

    # Menu
    for i in range(len(dists)):
        dist = dists[i]
        fh.write("<li><span style='cursor:pointer' onclick='show_div("+str(i+1)+")'>"+dist['title']+"</span>\n")

    fh.write("""
            </ul>

            <div style="margin-top:4em;text-align:center">
                <div style="margin-right:2em;display:inline;cursor:pointer;" onclick="show_next(-1)">
                <svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 1000 1000" width="2em" height="2em" enable-background="new 0 0 1000 1000" xml:space="preserve">
                <g><path d="M500,10.4C229.4,10.4,10,229.6,10,500c0,270.4,219.4,489.6,490,489.6S990,770.4,990,500C990,229.6,770.6,10.4,500,10.4z M638.9,526.9L404.4,761.2c-12,12-31.3,12-43.3,0c-12-12-12-31.3,0-43.3L579.2,500L361.1,282.1c-12-11.9-12-31.3,0-43.3c12-12,31.3-12,43.3,0l234.5,234.3c7.3,7.3,9.6,17.4,7.9,26.9C648.5,509.5,646.2,519.6,638.9,526.9z" transform="translate(1000 0) scale(-1 1)"/></g>
                </svg>
                </div>

                <div style="display:inline;cursor:pointer;" onclick="show_next(1)">
                <svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" viewBox="0 0 1000 1000" width="2em" height="2em" enable-background="new 0 0 1000 1000" xml:space="preserve">
                <g><path d="M500,10.4C229.4,10.4,10,229.6,10,500c0,270.4,219.4,489.6,490,489.6S990,770.4,990,500C990,229.6,770.6,10.4,500,10.4z M638.9,526.9L404.4,761.2c-12,12-31.3,12-43.3,0c-12-12-12-31.3,0-43.3L579.2,500L361.1,282.1c-12-11.9-12-31.3,0-43.3c12-12,31.3-12,43.3,0l234.5,234.3c7.3,7.3,9.6,17.4,7.9,26.9C648.5,509.5,646.2,519.6,638.9,526.9z"/></g>
                </svg>
                </div>
            </div>

            </td>
            <td style="background-color:#eee;vertical-align:top;padding:1em;">
        """)

    # Images
    img = None
    for i in range(len(dists)):
        dist = dists[i]
        img = dist['img']+'.png'
        css = 'display:none'
        fh.write("<div style='"+css+"' id='div"+str(i+1)+"'><img src='data:image/png;base64,")
        with open(dir+'/'+img, "rb") as image_file:
            fh.write(base64.b64encode(image_file.read()))
        fh.write("""'>
            </div>
            """)

    # Phantom image and summary numbers
    fh.write("""
            <div id="div0" style='display:none;'>
                <div id='placeholder' style='opacity:0'></div>
                <table class='summary'>
                    <tr><td>reads total <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_reads']) + """<td><td></tr>

                    <tr><td>reads processed <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_processed']) + """<td style="text-align:right;">""" 
                        + percent(dat['SN']['n_processed'],dat['SN']['n_reads']) + """<td></tr>

                    <tr><td>reads excluded <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_excluded']) + """<td style="text-align:right;">""" 
                        + percent(dat['SN']['n_excluded'],dat['SN']['n_reads']) + """<td></tr>

                    <tr><td class="indent">.. duplicate <td style="text-align:right;">""" 
                        + bignum(dat['SN']['n_flag_DUP']) + """<td style="text-align:right;">""" 
                        + percent(dat['SN']['n_flag_DUP'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. MQ&lt;"""+dat['LM']['n_min_mq']+"""  <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_min_mq']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_min_mq'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. soft-clipped bases &gt;"""+dat['LM']['n_max_soft_clips']+""" <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_max_soft_clips']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_max_soft_clips'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. secondary <td style="text-align:right;">""" 
                        + bignum(dat['SN']['n_flag_SECONDARY']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_flag_SECONDARY'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. supplementary <td style="text-align:right;">""" 
                        + bignum(dat['SN']['n_flag_SUPPLEMENTARY']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_flag_SUPPLEMENTARY'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. unmapped <td style="text-align:right;">""" 
                        + bignum(dat['SN']['n_flag_UNMAP']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_flag_UNMAP'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. QC failed <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_flag_QCFAIL']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_flag_QCFAIL'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. no barcode <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_no_BX']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_no_BX'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td class="indent">.. N in a barcode <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_barcode_seq_has_n']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_barcode_seq_has_n'],dat['SN']['n_reads']) +"""<td></tr>

                    <tr><td>barcodes <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_barcodes']) + """<td><td></tr>

                    <tr><td class="indent">.. number of good reads &lt;"""+dat['LM']['n_min_reads_per_barcode']+""" <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_min_reads_per_barcode']) + """<td style="text-align:right;">"""
                        + percent(dat['SN']['n_min_reads_per_barcode'],dat['SN']['n_barcodes']) +"""<td></tr>

                    <tr><td>read pruning<td style="text-align:right;"><td><td></tr>

                    <tr><td class="indent">.. read spacing &lt;"""+dat['LM']['n_min_read_spacing']+"""bp <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_min_read_spacing']) + """<td><td></tr>

                    <tr><td class="indent">.. read count in fragment &lt;"""+dat['LM']['n_min_reads_per_fragment']+""" <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_min_reads_per_fragment']) + """<td><td></tr>

                    <tr><td>fragments total <td style="text-align:right;">"""
                        + bignum(dat['SN']['n_fragments']) + """<td><td></tr>

                </table>
            </div>
            
            <div id="graph" style="position:relative"></div>

            </td></tr>
            </table>
            </div>
            <script type="text/javascript">
                var idx = 0;
                var max_idx = """ + str(len(dists)) + """
                function show_div(id)
                {
                    document.getElementById('graph').innerHTML = document.getElementById('div'+id).innerHTML; 
                    idx = id;
                }
                function show_next(delta)
                {
                    idx = idx + delta;
                    if ( idx < 0 ) idx = max_idx;
                    else if ( idx > max_idx ) idx = 0;
                    show_div(idx); 
                }
                document.getElementById('placeholder').innerHTML = document.getElementById("div1").innerHTML;
                show_div(0);
            </script>
            </body></html>
        """)
    fh.close()

#--------------------

plots = []
dat  = parse_file(fname)

for key in dists4:
    if key not in dat: continue
    if key not in graphs: graphs.append(key)
    plot_dist4(dat[key])

for key in dists:
    if key not in dat: continue
    if key not in graphs: graphs.append(key)
    plot_dist(dat[key])

for key in chr_maps:
    if key not in dat: continue
    if key not in graphs: graphs.append(key)
    plot_dist(dat[key])

for key in graphs:
    if key not in dat: continue
    plots.append(dat[key])

write_html(dir+'/'+"plots.html", plots)



