#COMMAND LINE ARGUMENT SHOULD BE DATE OF START OF OBSERVATIONS eg 20231213
scp DECamObserver@observer2.ctio.noao.edu:/user/DECamObserver/$1.qcinv ./2024A/$1/
git add ./2024A/$1/$1.qcinv
git commit -am 'adding qcinv'
git push 
