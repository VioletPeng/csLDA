#!/bin/bash

time ./src/lda -est -alpha 0.1 -beta 0.01 -tau 0.5 -epsilon 0.01 -ntopics 40 -niters 2000 -savestep 1000 -twords 20 -dfile ChineseData/Chinese.cs.test 
