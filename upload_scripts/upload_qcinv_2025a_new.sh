#!/bin/sh

#COMMAND LINE ARGUMENT SHOULD BE DATE OF START OF OBSERVATIONS YYMMDD eg 231213
scp  DECamObserver@observer4.ctio.noao.edu:/home/DECamObserver/20$1.qcinv  ./2025A/$1/
git add ./2025A/$1/20$1.qcinv
git commit -am 'adding qcinv'
git push 
