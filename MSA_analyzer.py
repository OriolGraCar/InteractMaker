
def analyse_MSA(sequence1,sequence2):
	count1 = 1
	count2 = 1
	resis1 = []
	resis2 = []
	for i in range(len(string1)):
		if string1[i] == '-':
			count1 -= 1
		if string2[i] == '-':
			count2 -= 1
		elif string1[i] != '-' and string2[i] != '-':
			resis1.append(i+count1)
			resis2.append(i+count2)
	return resis1,resis2



