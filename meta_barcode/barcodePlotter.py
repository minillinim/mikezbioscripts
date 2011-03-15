#!/usr/bin/python
import subprocess
import os
from optparse import OptionParser
import sys
from rpy2 import *
import rpy2.robjects as robjects
import re
from operator import itemgetter

###############################################################################
#
#    barcodePlotter.py
#
#    Cluster barcodes and produce heatmap images for barcode metrics
#    obtained using barcodeMers.pl
#    
#    Copyright (C) 2009 2010 Lauren Bragg, Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
# IMAGE Stuff
#
def openRImageDevice(label, r):
    r("png(filename='"+label+"', width=1000, height=700);")
    r('library(gplots);')
    

def closeRImageDevice(r):
    r("dev.off();")

def loadData(fileName, r):
    # print out the purdy picture
    r('windows <- read.table("'+fileName+'", header=TRUE, sep=",");')
    r('mymatrix <- as.matrix(windows[,2:dim(windows)[2]]);')
    r('data_columns <- dim(mymatrix)[2];')
    
def produceMap(reportFileName, num_clusts, breaks, commaBreak, isAverages, r):
    # print out the purdy picture
    
    # get the right breaks
    break_string = ''
    num_breaks = 0
    if(commaBreak is False):
        break_split = breaks.split(':')
        num_breaks = int(break_split[1]) 
        break_string = makeCommaBreak(int(break_split[0]),(int(break_split[0]) / int(break_split[1])))
    else:
        break_string = breaks
        breaks_split = breaks.split(',')
        num_breaks = len(breaks_split) - 1
    
    # make a palette
    r('rgbs_G <- cbind(seq(0,180, length=256), seq(0,255, length=256), seq(0,180, length=256)) / 255;')
    r('cols_G <- rgb(rgbs_G[,1], rgbs_G[,2], rgbs_G[,3]);')
    col_step = int((255/num_breaks))
    
    # get the first colour
    r('fred <-cols_G[1];')
    bob = r['fred']
    col_string='"'+str(bob[0])+'"'
    cur_col = col_step;
    
    # now get the rest
    for i in range (1,num_breaks):
        r('fred <-cols_G['+str(cur_col)+'];')
        cur_col+=col_step
        bob = r['fred']
        col_string=col_string+',"'+str(bob[0])+'"'
    
    # fix the side bar colours
    avail_cols = ['"#BBBBBB"','"#333333"']
    side_cols_string = 'rep('+avail_cols[1]+', numClusters[1])' 
    
    for i in range(2, int(num_clusts)+1):
        side_cols_string = side_cols_string+',rep('+avail_cols[i%2]+', numClusters['+str(i)+'])'
    
    # load colors into R
    r('colors <- c('+col_string+');')
    r('breaks <- c('+break_string+')')
    r('rsc <- c('+side_cols_string+')')
    
    # draw the heatmap
    r('heatmap.2(sorted_comb[,(-1*max_col)], scale="none", col=colors, breaks=breaks, trace="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, labCol=FALSE, labRow=FALSE, rowsep=c(705,2384), revC=FALSE, RowSideColors=rsc);')

    if(isAverages is False):
        # make a report formulti contig windows
        r('counts <- c()')
        r('labels <- c()')
        r('cumu_string <- c()')

        #Get all the unique rownames (minus cut number on the end)
        r('rns <- rownames(combined);')
        bob=r['rns']
        for i in range(0,len(bob)):
            bob[i] = re.sub(r'__.*', '',bob[i])
        
        r('contig_names <- unique(rownames(combined))')
        r('data_rows <- length(contig_names);')
        r('contig_exp <- rownames(combined)')
        fred = r['contig_names']
        billy = r['contig_exp']
        
        l_zzz = ['nothing']
        l_123 = []
        for i in range(1,int(num_clusts)+1):
            l_zzz.append('0')
            l_123.append(str(i))
        c_zzz = robjects.StrVector(l_zzz)
        c_123 = robjects.IntVector(l_123)
        r('contig_calls <- '+c_zzz.r_repr())
        r('for(i in 1:length(contig_names)) \
        { \
            nameI <- contig_names[i]; \
            myRow <- table(c(combined[which(contig_exp == nameI),max_col],'+c_123.r_repr()+')); \
            fc <- formatC(c(nameI, myRow), format="d"); \
            cumu_string <- cbind(cumu_string, fc); \
        }')
        out_p = r['cumu_string']
        num_rows = r['data_rows']
        
        # print out a quality report
        report_file = open(reportFileName, 'wb')
        for i in range(0, num_rows[0]):
            for j in range(0, int(num_clusts)+1):
                no_ws = re.sub(r'\s', '', out_p[i*(int(num_clusts)+1) + j])
                report_file.write(no_ws+",")
            report_file.write("\n")
        report_file.close()
    else:
        # make report for contig averages
        r('clust_names = rownames(combined)')
        r('clust_values = combined[,max_col]')
        cv = r['clust_values']
        rn = r['clust_names']
        bobby = dict()
        for i in range(0,len(cv)):
            bobby[rn[i]] = cv[i]
        
        # now bobby is a key value dictionary, so sort it
        bobby_sorted = sorted(bobby.items(), key=itemgetter(1))
        
        # output
        report_file = open(reportFileName, 'wb')
        curr_clust = 1;
        report_file.write("*************************************\n\tContigs in cluster: "+str(curr_clust)+"\n*************************************\n\n")
        
        for i in  range (0, len(bobby_sorted)):
            if(curr_clust != bobby_sorted[i][1]):
                # new cluser
                curr_clust = bobby_sorted[i][1]
                report_file.write("\n*************************************\n\tContigs in cluster: "+str(curr_clust)+"\n*************************************\n\n")
            
            report_file.write(bobby_sorted[i][0]+"\n")
        
        report_file.close()
            
        
#
# Cluster Testing
#
def clustRTest(cent_range, r):
    #using clValid
    
    r('library(clValid);')
    r('intern <- clValid(mymatrix, '+cent_range+', clMethods = c("hierarchical", "kmeans"), validation = "internal");')

    # test 1
    r('pdf("internal_validation_clusters.pdf");')
    r('op <- par(no.readonly = TRUE);')
    r('par(mfrow = c(2, 2), mar = c(4, 4, 3, 1));')
    r('plot(intern, legend = FALSE);')
    r('plot(nClusters(intern), measures(intern, "Dunn")[, , 1], type = "n", axes = F, xlab = "", ylab = "");')
    r('legend("center", clusterMethods(intern), col = 1:9, lty = 1:9,pch = paste(1:9));')
    r('par(op);')
    r('dev.off();')

    # test 2
    r('pdf("stability_validation_clusters.pdf");')
    r('stab <- clValid(mymatrix, '+cent_range+', clMethods = c("hierarchical", "kmeans"), validation = "stability");')
    r('optimalScores(stab);')
    r('par(mfrow = c(2, 2), mar = c(4, 4, 3, 1));')
    r('plot(stab, measure = c("APN", "AD", "ADM"), legend = FALSE);')
    r('plot(nClusters(stab), measures(stab, "APN")[, , 1], type = "n", axes = F, xlab = "", ylab = "");')
    r('legend("center", clusterMethods(stab), col = 1:9, lty = 1:9,pch = paste(1:9));')
    r('par(op);')
    r('dev.off();')
    
    # now get a tree for H-clustering
    r('pdf("heirachical_clusters.pdf")')
    r('dmat <- dist(mymatrix);')
    r('mytree <- hclust(dmat);')
    r('plot(mytree);')
    r('par(op);')
    r('dev.off();')
    

#
# Clustering K means
#
def clustKMeans(centers, r):
    # glorified sorting algorithm
    r('rownames(mymatrix) <- windows[,1];')
    r('run <- kmeans(mymatrix, centers='+centers+', iter.max=20);')
    r('combined <- cbind(mymatrix, run$cluster);')
    r('numClusters <- table(run$cluster);')
    r('max_col <- data_columns+1;')
    r('sorted_comb <- combined[order(combined[,max_col]), 1:data_columns];')

def clustHClust(centers, r):
    # glorified sorting algorithm
    r('rownames(mymatrix) <- windows[,1];')
    r('dmat <- dist(mymatrix);')
    r('mytree <- hclust(dmat);')
    r('clusters <- cutree(mytree, k='+centers+');')    
    r('combined <- cbind(mymatrix, clusters);')
    r('numClusters <- table(clusters);')
    r('max_col <- data_columns+1;')
    r('sorted_comb <- combined[order(combined[,max_col]), 1:data_columns];')

def makeCommaBreak(maxb, sepb):
    # make a comma break pattern from max and sep parameters
    break_string = '0'
    current_break = sepb;
    while (current_break <= maxb):
        break_string = break_string+','+ str(current_break)
        current_break += sepb
    return break_string
    

def printAtStart():
    print "---------------------------------------------------------------- "
    print "barcodePlotter.py"
    print "Copyright (C) 2009, 2010 Lauren bragg, Michael Imelfort\n"
        
    print "This program comes with ABSOLUTELY NO WARRANTY;"
    print "This is free software, and you are welcome to redistribute it"
    print "under certain conditions: See the source for more details."
    print "---------------------------------------------------------------- "

#
# Entry sub. Parse vars and call parseSamBam
#
if __name__ == '__main__':

    # intialise the options parser
    parser = OptionParser("\n\n %prog -b barcodes -c { num_clusters | cluster_range_start:cluster_range_end } -l image label [-H] [-a] [-S]")
    parser.add_option("-b", "--barcode_fileName", type="string", dest="barcodeFileName", help="Specify a name for the barcodes CSV file")
    parser.add_option("-c", "--clusters", type="string", dest="numClusts", help="The number of different clusters you expect to have. OR the range you wish to test (for k-means only) a:b  (a > 1)")
    parser.add_option("-r", "--report", type="string", dest="reportName", help="The file name for the report to be printed to [default: barcode_plotz.log]")
    parser.add_option("-H", "--hierarchical", action="store_true", dest="doHClust", help="Set this to use hierarchical clustering [default: k-means]")
    parser.add_option("-B", "--breaks", type="string", dest="imageBreaks", help="Specify the breaks for colouring '0,b1[,bn]' OR 'max:number_of_breaks' will be segregated linearly ( [default: 16:8]")
    parser.add_option("-l", "--image_name", type="string", dest="imageName", help="Specify a label for the image")
    parser.add_option("-a", "--averages", action="store_true", dest="isAverages", help="Set this if the barcode file contains whole contigs window averages [default: false]")
    parser.add_option("-S", "--Silent", action="store_true", dest="doSilent", help="Print nothing extra! [default: false]")
    
    # get and check options
    (opts, args) = parser.parse_args()
    
    if (opts.doSilent is None):
        printAtStart()

    if (opts.barcodeFileName is None):
        print ('\nERROR: You need to specify a .csv file to parse\n')
        parser.print_help()
        sys.exit(1)
        
    # now check to see if we are testing a range or clustering nicely?
    if (opts.numClusts is None):
        print ('\nERROR: You need to specify the number of clusters\n')
        parser.print_help()
        sys.exit(1)
    if(re.search(':', opts.numClusts) is None):
        # straight up custer and print
        only_clust = False
    else:
        only_clust = True
        
    if((opts.imageName is None) and (only_clust is False)):
        print ('\nERROR: You need to specify an image name\n')
        parser.print_help()
        sys.exit(1)
        
    if(opts.isAverages is None):
        doAves = False
    else:
        doAves = True

    # get a log file happenin'
    barcode_log = "barcode_plotz.log"
    if(opts.reportName):
        barcode_log = opts.reportName

    #check that the breaks are kosher
    breaks = "16:8"
    comma_break = False
    if(opts.imageBreaks):
        if(re.search(':', opts.imageBreaks)):
            # the user has specified a max/min break
            test_mod_split = opts.imageBreaks.split(':')
            if(0 != (int(test_mod_split[0]) % int(test_mod_split[1]))):
                print ('\nERROR: The number of breaks must evenly divide the max count!\n')
                sys.exit(1)  
            breaks=opts.imageBreaks
        else: 
            #the user has input the comma form of breaks.
            comma_break = True
            breaks = opts.imageBreaks
    
    # open an R instance and load the data
    r = robjects.r
    loadData(opts.barcodeFileName, r)
    
    if(only_clust):
        # do clustering only!
        print "\n\n***Cluster testing.\nyou will need to look at the following files:\n * heirachical_clusters.pdf\n * internal_validation_clusters.pdf\n * stability_validation_clusters.pdf\n\n"
        clustRTest(opts.numClusts, r)
        sys.exit(0)

    # clustering
    if(opts.doHClust is None):
        # kmeans
        clustKMeans(opts.numClusts, r)
    else:
        clustHClust(opts.numClusts, r)

    # make the heatmap
    openRImageDevice(opts.imageName, r)
    produceMap(barcode_log, opts.numClusts, breaks, comma_break, doAves, r)
    closeRImageDevice(r)
