#!/usr/bin/python
import subprocess
import os
from optparse import OptionParser
import sys
from rpy2 import *
import rpy2.robjects as robjects
###############################################################################
#
#    pipelinePlotter.py
#    Produce images plots from the parsed csv made in sam2PairPlotCSV
#    Copyright (C) 2010 Adam Skarshewski, Michael Imelfort
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
# Directionality
###############################################################################
#  0  DD  <---  --->   if the first read disagrees with ref then the second is ahead and in the reverse orientation (wrt the first)
#  3  DA  <---  <---   if the first read disagrees with ref then the second is ahead and in the same orientation
#  2  AD  --->  <---   if the first read agrees with ref then the second is ahead and in the reverse orientation
#  3  AA  --->  --->   if the first read agrees with ref then the second is ahead and in the same orientation
#  4  Single hit: this read agrees with it's ref
#  5  Single hit: this read disagrees with it's ref
#  6  Different references: this read agrees with it's ref
#  7  Different references: this read disagrees with it's ref
###############################################################################

#
# Stuff
#
def openRImageDevice(label, r):
    png_filename = label + ".png"
    r("png(filename='"+png_filename+"', width=1000, height=700);")
    r("par( mfrow= c(4, 1));")


def closeRImageDevice(r):
    r("dev.off();")

def plotPairs(fileName, r):

    # set the transform to a log
    r("transformation = function (X) { log(X) }")
    #r("transformation = function (X) { X }")

    r("par(mar=c(1, 4, 0, 2) + 0.1);")

    # readin and transform the file
    r("x = read.csv(\""+fileName+"\")")
    print "{{{" + str(r("maxsert = max(x$insertsize) * 1.5")) + "}}}"
    #r("scaling = (10^floor(log(max(x$position),10)-1))")
    #r("scaffoldsize = max(x$position)")
    #r("roundedscaffoldsize = scaling*(floor(max(x$position)/scaling)+1)")

    # check we don't have any zeros which will kill the log transform
    maxsert = r['maxsert']
    print "Max insert is: {{" + str(maxsert) + "}}"
    if(maxsert[0] == 0):
            print "Max insert is 0"
            sys.exit(1)
    r("ylims=c(-5,transformation(maxsert))")

#################
# Data Containers
#################
    r("opposed = list(position = x$position[x$direction == 0], insertsize = x$insertsize[x$direction == 0])")
    r("facing = list(position = x$position[x$direction == 2], insertsize = x$insertsize[x$direction == 2])")
    r("same = list(position = x$position[x$direction == 3], insertsize = x$insertsize[x$direction == 3])")
    r("forwards = list(position = x$position[x$direction == 4], insertsize = x$insertsize[x$direction == 4])")
    r("reverses = list(position = x$position[x$direction == 5], insertsize = x$insertsize[x$direction == 5])")
    r("forwardjumpers = list(position = x$position[x$direction == 6], insertsize = x$insertsize[x$direction == 6])")
    r("reversejumpers = list(position = x$position[x$direction == 7], insertsize = x$insertsize[x$direction == 7])")

    # find the mode of the largest occurring insert type
    r("op_len = length(opposed$position);")
    r("fa_len = length(facing$position);")
    r("sa_len = length(same$position);")
    max = r['op_len']
    max_str = "opposed"
    fa = r['fa_len']
    sa = r['sa_len']
    if(fa[0] > max[0]):
        max[0] = fa[0]
        max_str = "facing"
    if(sa[0] > max[0]):
        max_str = "same"
    r("ins_mode = mean("+max_str+"$insertsize)")

#################
# PowerPlot
#################
    r("plot(c(0,roundedscaffoldsize), c(0,transformation(maxsert)), xlim=c(0,roundedscaffoldsize), ylim=ylims, pch='*', xlab='', ylab='(ln) Insert Size', axes='FALSE', col='white'); \
       axis(2, c(0, floor(transformation(ins_mode))), labels=c('0',as.character(floor(ins_mode))), las='2'); \
       axis(1, c(0, roundedscaffoldsize/5, 2*roundedscaffoldsize/5, 3*roundedscaffoldsize/5, 4*roundedscaffoldsize/5, roundedscaffoldsize)); \
       grid(); \
       points(opposed$position, transformation(opposed$insertsize), pch='.', col='red'); \
       points(facing$position, transformation(facing$insertsize), pch='.', col='black'); \
       points(same$position, transformation(same$insertsize), pch='.', col='blue'); \
       points(forwards$position, transformation(forwards$insertsize) - 1, pch='>', col='purple'); \
       points(reverses$position, transformation(reverses$insertsize) - 2, pch='<', col='orange'); \
       points(forwardjumpers$position, transformation(forwardjumpers$insertsize) - 3, pch='>', col='green'); \
       points(reversejumpers$position, transformation(reversejumpers$insertsize) - 4, pch='<', col='magenta');")


#
# insert size histogram
#
    r("par(mar=c(2, 4, 1, 2) + 0.1);")
    r("mean_mult = 2")
    r("phil <- x$insertsize[which(x$insertsize < mean_mult*ins_mode)];")
    r("phil <- phil[which(phil > 1)];")
    r("hist(phil, breaks=200, col='red', xlab='', ylab='Insert Size freq.', main='', xlim=c(0,mean_mult*ins_mode));")
#    r("hist(c(opposed$insertsize[opposed$insertsize], facing$insertsize[facing$insertsize]), breaks=100, col='red', xlab='', ylab='Insert Size freq.', main='', xlim=c(0,1.5*ins_mode)); \
#       histobj = hist(facing$insertsize[facing$insertsize < 1000], breaks=100, plot = FALSE); \
#       lines(histobj, col='black');")

#
# plot coverage
#
def plotCoverage(fileName, modeFileName, r):

    r("par(mar=c(1, 4, 0, 2) + 0.1);")

    r("cov_data = read.csv(\""+fileName+"\")")
    r("max_cov = max(cov_data$coverage) + 5")
    r("scaling = (10^floor(log(max(cov_data$position),10)-1))")
    r("scaffoldsize = max(cov_data$position)")
    r("roundedscaffoldsize = scaling*(floor(max(cov_data$position)/scaling)+1)")

    r("mode_cov_data = read.csv(\""+modeFileName+"\");")
    r("mode_cov = max(mode_cov_data$mode_coverage);")
    r("mean_cov = mean(cov_data$coverage)")

    # check we don't have any zeros which will kill the log transform
    max_cov = r['max_cov']
    if(max_cov[0] == 0):
        print "Max coverage is 0!"
        sys.exit(1)

    r("ylims=c(0,max_cov);")
    r("plot(c(0,roundedscaffoldsize), c(0,max_cov), xlim=c(0,roundedscaffoldsize), ylim=ylims, pch='*', xlab='', ylab='Vertical Coverage', axes='FALSE', yaxt='n', col='white'); \
       axis(1, c(0, roundedscaffoldsize/5, 2*roundedscaffoldsize/5, 3*roundedscaffoldsize/5, 4*roundedscaffoldsize/5, roundedscaffoldsize)); \
       axis(2, c(0, floor(max_cov/3), floor(2*max_cov/3), floor(max_cov))); \
       grid(); \
       lines(cov_data$position, cov_data$coverage, pch='.', col='#FF034A'); \
       lines(c(0,roundedscaffoldsize), c(mode_cov, mode_cov), col='black'); \
       lines(c(0,roundedscaffoldsize), c(mean_cov, mean_cov), col='green')")

def plotKmerHeatMap(fileName, r):

    # first we need to find out how many colums we have in the
    # file and then use these ar labels
    csv_file = open(fileName, 'r')
    line = csv_file.readline()
    line = line.replace('"','')
    line = line.rstrip('\n')
    cols = line.rsplit(',')
    csv_file.close()
    num_cols = len(cols);

    # we need this so we can see the bottom axis ticks
    r("par(mar=c(4, 4, 0, 2) + 0.1);")

    # read in and transform the file
    #r("transformation = function (X) { log(X) }")
    r("transformation = function (X) { X }")
    r("hm_data = read.csv(\""+fileName+"\")")
    #r("scaling = (10^floor(log(max(hm_data$position),10)-1))")
    #r("scaffoldsize = max(hm_data$position)")
    #r("roundedscaffoldsize = scaling*(floor(max(hm_data$position)/scaling)+1)")

    # get nice green-red and red-fade pallettes
    r("rgbs_R <- cbind(255, seq(255,0, length=256), 0) / 255; \
       rgbs_G <- cbind(seq(0,255, length=256), 255, 0) / 255; \
       rgbs_RG = (rgbs_R * rgbs_G ); \
       cols_RG <- rgb(rgbs_RG[,1], rgbs_RG[,2], rgbs_RG[,3]); \
       rgbs_R <- cbind(255, seq(255,32, length=256), 0) / 255; \
       cols_R <- rgb(rgbs_R[,1], rgbs_R[,2], rgbs_R[,3]);")

    r("plot(c(0,roundedscaffoldsize), c(0,"+str(num_cols)+"), xlim=c(0,roundedscaffoldsize), ylim=c(0,"+str(num_cols)+"), pch='*', xlab='', ylab='Kmer abundance', axes='FALSE', col='white');")
    r("grid();")
    #return

    for i in range (1, num_cols - 1):
        r("total_len = length(hm_data$position)")
        r("bob = hm_data$position-hm_data$position+"+str(i)+"")
        # find the y limits this will be the maximum of the fist column...
        r("max_y = transformation(max(hm_data$"+cols[i]+") + 10)")
        r("col_scaler=255/max_y")
        r("scaled_col_ids = hm_data$"+cols[i]+" * col_scaler")
        r("plot_cols = cols_R[scaled_col_ids]")
        r("points(hm_data$position, bob, pch='|', col=rep(plot_cols, each=1, length.out=length(hm_data$position)))")

    r("col_scaler=255/5")
    r("scaled_col_ids = hm_data$"+cols[num_cols-1]+" * col_scaler")
    r("plot_cols = cols_RG[scaled_col_ids]")
    r("points(hm_data$position, bob+1, pch='|', col=rep(plot_cols, each=1, length.out=length(hm_data$position)))")

    # draw the purdy axis
    r("axis(1, c(0, roundedscaffoldsize/5, 2*roundedscaffoldsize/5, 3*roundedscaffoldsize/5, 4*roundedscaffoldsize/5, roundedscaffoldsize))")
    r("axis(2, c(1, 2, 3, 4,5), labels=c('3','4','8','12','GC'), las='2');")

def plotGC(fileName, r):

    # first we need to find out how many colums we have in the
    # file and then use these ar labels
    csv_file = open(fileName, 'r')
    line = csv_file.readline()
    line = line.replace('"','')
    line = line.rstrip('\n')
    cols = line.rsplit(',')
    csv_file.close()
    num_cols = len(cols);

    # we need this so we can see the bottom axis ticks
    r("par(mar=c(4, 4, 0, 2) + 0.1);")

    # read in and transform the file
    #r("transformation = function (X) { log(X) }")
    r("transformation = function (X) { X }")
    r("hm_data = read.csv(\""+fileName+"\")")
    #r("scaling = (10^floor(log(max(hm_data$position),10)-1))")
    #r("scaffoldsize = max(hm_data$position)")
    #r("roundedscaffoldsize = scaling*(floor(max(hm_data$position)/scaling)+1)")

    # get nice green-red and red-fade pallettes
    r("rgbs_R <- cbind(255, seq(255,0, length=256), 0) / 255; \
       rgbs_G <- cbind(seq(0,255, length=256), 255, 0) / 255; \
       rgbs_RG = (rgbs_R * rgbs_G ); \
       cols_RG <- rgb(rgbs_RG[,1], rgbs_RG[,2], rgbs_RG[,3]); \
       rgbs_R <- cbind(255, seq(255,32, length=256), 0) / 255; \
       cols_R <- rgb(rgbs_R[,1], rgbs_R[,2], rgbs_R[,3]);")

    r("plot(c(0,roundedscaffoldsize), c(0,"+str(num_cols)+"), xlim=c(0,roundedscaffoldsize), ylim=c(0,"+str(num_cols)+"), pch='*', xlab='', ylab='Kmer abundance', axes='FALSE', col='white');")
    r("grid();")
    #return

    for i in range (1, num_cols - 1):
        r("total_len = length(hm_data$position)")
        r("bob = hm_data$position-hm_data$position+"+str(i)+"")
        # find the y limits this will be the maximum of the fist column...
        r("max_y = transformation(max(hm_data$"+cols[i]+") + 10)")
        r("col_scaler=255/max_y")
        r("scaled_col_ids = hm_data$"+cols[i]+" * col_scaler")
        r("plot_cols = cols_R[scaled_col_ids]")
        r("points(hm_data$position, bob, pch='|', col=rep(plot_cols, each=1, length.out=length(hm_data$position)))")

    r("col_scaler=255/5")
    r("scaled_col_ids = hm_data$"+cols[num_cols-1]+" * col_scaler")
    r("plot_cols = cols_RG[scaled_col_ids]")
    r("points(hm_data$position, bob+1, pch='|', col=rep(plot_cols, each=1, length.out=length(hm_data$position)))")

    # draw the purdy axis
    r("axis(1, c(0, roundedscaffoldsize/5, 2*roundedscaffoldsize/5, 3*roundedscaffoldsize/5, 4*roundedscaffoldsize/5, roundedscaffoldsize))")
    r("axis(2, c(1, 2, 3, 4,5), labels=c('3','4','8','12','GC'), las='2');")

#
# Entry sub. Parse vars and call parseSamBam
#
if __name__ == '__main__':

    # intialise the options parser
    parser = OptionParser("\n\n %prog -c csvFileName [-l image label]")
    parser.add_option("-p", "--pair_fileName", type="string", dest="pairCSVFileName", help="Specify a name for the pair plot CSV file")
    parser.add_option("-c", "--cov_fileName", type="string", dest="covCSVFileName", help="Specify a name for the coverage CSV file")
    parser.add_option("-m", "--mode_cov_fileName", type="string", dest="modeCovCSVFileName", help="Specify a name for the mode coverage CSV file")
    parser.add_option("-H", "--heatmap_fileName", type="string", dest="heatCSVFileName", help="Specify a name for the heat map CSV file")
    parser.add_option("-l", "--image_label", type="string", dest="imageLabel", help="Specify a label for the image")

    # get and check options
    (opts, args) = parser.parse_args()
    if (opts.pairCSVFileName is None):
        print ('You need to specify a .csv file to parse for paired plotting')
        parser.print_help()
        sys.exit(1)
    if (opts.covCSVFileName is None):
        print ('You need to specify a .csv file to parse for coverage plotting')
        parser.print_help()
        sys.exit(1)
    if (opts.modeCovCSVFileName is None):
        print ('You need to specify a .csv file to parse for mode coverage plotting')
        parser.print_help()
        sys.exit(1)
    if (opts.heatCSVFileName is None):
        print ('You need to specify a .csv file to parse for heatmap plotting')
        parser.print_help()
        sys.exit(1)

    # override defaults
    if(opts.imageLabel is None):
        label = 'ADD_LABEL'
    else:
        label =  opts.imageLabel

    # do stuff
    r = robjects.r

    openRImageDevice(label, r)
    plotCoverage(opts.covCSVFileName, opts.modeCovCSVFileName, r)
    plotKmerHeatMap(opts.heatCSVFileName,  r)
    plotPairs(opts.pairCSVFileName, r)
    closeRImageDevice(r)
