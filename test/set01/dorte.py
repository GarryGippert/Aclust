#! /usr/bin/python3
''' test some combinations of input '''
import os
def doc(cmd):
	s = os.system(cmd)
	if s > 0:
		raise Exception(s, cmd)
	print("OK", s, cmd)

cmd = "../../src/bclust -p one -alf AA.aln AB.aln BB.aln > one.out 2> one.err"
doc(cmd)
print ("bclust runs with several ALF input files")

cmd = "diff one.dree.txt expect.dree.txt"
doc(cmd)
print ("bclust produces expected output with several ALF input files")

cmd = "cat AA.aln AB.aln BB.aln > AAABBB.aln"
doc(cmd)

cmd = "../../src/bclust -p two -alf AAABBB.aln > two.out 2> two.err"
doc(cmd)
print ("bclust runs with concatenated ALF input")

cmd = 'diff two.dree.txt expect.dree.txt'
doc(cmd)
print ("bclust produces expected output with concatenated ALF input")
