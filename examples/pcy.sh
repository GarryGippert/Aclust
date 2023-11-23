echo start `date`
echo "doing ../bin/aclust -s ../dat/BLOSUM62.dat pcy.fa > pcy.aclust.out 2> pcy.aclust.err"
../bin/aclust -s ../dat/BLOSUM62.dat pcy.fa > pcy.aclust.out 2> pcy.aclust.err
echo done `date`
