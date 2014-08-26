#!/bin/sh

# sendtotext.sh
# 
#
# Created by Tim Shokair on 3/22/11.
# Copyright 2011 University of Pennsylvania. All rights reserved.

#FILES =/Users/timshokair/Desktop/analysis*
for f in *.root
do 
	echo $PWD$f >>/Users/timshokair/Desktop/analysis/filenames.txt
	
done
