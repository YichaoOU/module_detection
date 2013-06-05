import sys
import os
import numpy as np
import scipy as sp
from decimal import Decimal
import re
'''
Global settings: PatScan path, Threshold for minimum number of hits, zScore threshold, largest Distance between two motifs

Class:
	statistics: 	calculate if the 'occurence of certain distance' (denote as x) is significant
		__init__(self)
		sigP(self,n,x)	calculate the p-value based on binomial distribution
		choose(self,n,k)
		oOverE(self,n,x)
		zScore(self,var,mean,x)
	
	motifs:			Read the motif file, create Modules based on the given motifs
		__init__(self,motif_file)
		readFile(self)
		createModules(self,num)
			num is the number of motifs used to create module, e.g. if num = 3, then it will use all possible combination of 3 motifs to constrcut the module
		displayMotif(self)	print out to screen
		displayModule(self)
		setModule(self)
			contruct the module based on the order of motifs (top-down) in the motif file
			this method will only create one module
	
	modules:		call PatScan to search patterns
		__init__(self,module_array)	module_array is the motifs.module from class 'motifs'
		searchModule(self)			call PatScan to search patterns
		testStatistics(self)
		display(self,thresHits)		thresHits is a threshold for the minimun number of hits of a given module

motif file:
	one motif(just A,C,G,T) occupy one line, example:
		ATGTGG
		AATGTTGGTC
		ATGC

typical use:
	import this script
	use motifs to read your motif file
	use modules to display significant module
	
	1. de novo searching (When you have lots of motifs and you have no idea what the module will looks like)
	
		test = motifs(motif_file)
		test.readFile()
		test.createModules()
	
		test_module = modules(test.module)
		test_module.searchModule()
		test_module.testStatistics()
		test_module.display(20)
		
	2. target module searching (When you have pre-knowledge. e.g. you know the module can only be A->B->N, not N->B->A)
		
		for each motif file:
		
			test = motifs(motif_file)
			test.readFile()
			test.setModule() ----- This is the only difference
	
			test_module = modules(test.module)
			test_module.searchModule()
			test_module.testStatistics()
			test_module.display(20)
		
'''
PatScan = "./scan_for_matches -c " 
zThres = 3.5 # a powerful default value
largestDistance = 60
distant = " 1...60 " # distant between each motif (largestDistance)
thresHits = 20

def nextLeaf(array,Length,k):
	for i in range(Length):
		if array[i]<k:
			array[i] += 1
			return array
		array[i] = 0

def uniqueCount(array):
	output =[]
	for item in range(1,largestDistance+1):
		output.append(array.count(item))
	return output

	
	

class statistics:
	
	def __init__(self):
		self.pValue = 1.0
		self.prob = 1./60. #This is depend on your distant setting, if you want to calculate p-value, this is needed
		self.zValue = 0.0
		
	def sigP(self,n,x): #This is a test for calculate the P-value, which is not very fast
		if x == 0:
			self.pValue = 1.00
		else:
			for i in range(x,n+1):
				a=Decimal(self.prob**i)
				b=Decimal((1-self.prob)**(n-i))
				self.pValue = self.pValue + Decimal(self.choose(n,i))*a*b
				
	def choose(self,n, k): # Mathmatical funtion
		if 0 <= k <= n:
			ntok = 1
			ktok = 1
			for t in xrange(1, min(k, n - k) + 1):
				ntok *= n
				ktok *= t
				n -= 1
			return ntok // ktok
		else:
			return 0
			
	def oOverE(self,n,x): # O/E score has bias
		if x > 4*n*self.prob:
			self.pValue = 0.0
	
	def zScore(self,var,mean,x): # zScore is fast and powerful
		if var == 0:
			self.zValue = 0.0
		else:
			self.zValue = abs(mean-x)/np.sqrt(var)
	
class motifs:
	motif = []
	module = []
	name = []
	moduleName = []

	def __init__(self,motif_file):
		self.file = motif_file
		self.distant=distant
		
	def createModules(self,num):
		array = num*[0]
		while True:
			module_current = self.motif[array[0]]
			for i in range(1,num):
				module_current += self.distant+self.motif[array[i]]
			self.module.append(module_current)
			array = nextLeaf(array,num,len(self.motif)-1)
			if array.count(len(self.motif)-1) == num:
				break
		module_current = self.motif[array[0]]
		for i in range(1,num):
			module_current += self.distant+self.motif[array[i]]
		self.module.append(module_current)
		
	def assemblyModules(self,module_array):
		new_module = []
		for module in module_array:
			for motif in self.motif:
				newModule = module + distant + motif
				new_module.append(newModule)
		return new_module
		
	def readFile(self):
		constraint = "[1,0,0]" # allow 1 mismatch
		for line in open(self.file).readlines():
			line=line.split()
			for item in line:
				self.motif.append(item+constraint)
	
	def readPWM(self,thres):
		line = open(self.file).readlines()
		
		for i in range(0,len(line)-4,6):
			currentLine = line[i].split()
			self.name.append(currentLine[0])
			array = []
			max_total_Score = 0.0
			for j in range(1,5):
				col=[]
				currentLine = line[i+j].split()
				max_single_Score = 0.0
				for x in range(1,len(currentLine)):
					temp = float(currentLine[x])*100
					if temp>max_single_Score:
						max_single_Score = temp
					col.append(str(int(temp)))
				max_total_Score += max_single_Score
				array.append(col)
			constraint = " >"+str(int(thres*max_total_Score))+" "
			matrix = np.array(array)
			[m,n]=matrix.shape
			pwmPattern=[]
			for j in range(n):
				pos = ""
				for x in range(3):
					pos += matrix[x,j]+","
				pos = "("+pos+matrix[3,j]+")"
				pwmPattern.append(pos)
			motifPattern="{"
			for item in range(len(pwmPattern)-1):
				motifPattern+=pwmPattern[item]+","
			motifPattern+=pwmPattern[len(pwmPattern)-1]+"}" + constraint
			self.motif.append(motifPattern)
	
	def displayMotif(self):
		for i in range(len(self.motif)):
			print self.motif[i]
	
	def displayModule(self):
		for i in range(len(self.module)):
			print self.module[i]
			
	def setModule(self):
		constraint = "[1,0,0]" # allow 1 mismatch
		module_current = self.motif[0]+constraint
		for i in range(1,len(self.motif)):
			module_current += self.distant+self.motif[i]+constraint			
		self.module.append(module_current)
	
	def displayName(self):
		for item in self.moduleName:
			print item

class modules:

	module = []	
	sigRes = []
	hits   = []
	num_dis = 0
	
	def __init__(self,module):
		self.module = module
		self.result = []
		self.dis_result = []
		
	def searchModule(self,Promoter):
		self.Promoter = Promoter
		for item in self.module:
			f = open('pattern_file','w')
			f.write(item)
			f.close()
			command = PatScan+"pattern_file < " + Promoter
			self.result.append(os.popen(command).readlines())

	def testStatistics(self):
		for item in range(len(self.result)):
			line = self.result[item]
			self.hits.append(len(line))
			if line == []:
				self.sigRes.append("null")
				continue
			array = []
			for i in range(len(line)):
				if i%2 != 0:
					tempLine = line[i].split()
					col=[]
					for j in range(len(tempLine)):
						if j%2 != 0:
							col.append(len(tempLine[j]))
					array.append(col)
					self.num_dis = len(col)
			matrix = np.array(array)
			aboveP = []
			for i in range(self.num_dis):
				list = matrix[:,i].tolist()
				uList = uniqueCount(list)
				variance = np.var(uList)
				mean     = np.mean(uList)
				sigDis = []
				for j in range(1,largestDistance+1):
					a = statistics()
					a.zScore(variance,mean,list.count(j))
					z_val = a.zValue
					if z_val >= zThres:
						sigDis.append(j)
				aboveP.append(sigDis)
			self.sigRes.append(aboveP)
	
	def testZScore(self):
		output_module=[]
		for item in range(len(self.result)):
			line = self.result[item]
			if line == []:
				continue
			array = []
			for i in range(len(line)):
				if i%2 != 0:
					tempLine = line[i].split() # it looks like: AATGTGT A ATGTGC
					j = len(tempLine) - 2 # The last but one is the dis, and you need to get its length
					array.append(len(tempLine[j]))
			uList = uniqueCount(array)
			variance = np.var(uList)
			mean     = np.mean(uList)
			for j in range(1,largestDistance+1):
				a = statistics()
				a.zScore(variance,mean,array.count(j))
				z_val = a.zValue
				if z_val >= zThres:
					Original_module = self.module[item]
					pre_dis = 0
					post_dis = j+1
					if (j <= 1):
						pre_dis = 1
					else:
						pre_dis = j-1
					new_dis = " "+str(pre_dis)+"..."+str(post_dis)+" "
					new_sig_module = re.sub(distant,new_dis,Original_module)
					output_module.append(new_sig_module)
		return output_module
	
	def display(self,thresHits):
		for i in range(len(self.module)):
			flag = 0
			for j in range(len(self.sigRes[i])):
				if self.sigRes[i] == "null" or self.hits[i] < thresHits or len(self.sigRes[i][j]) == 0:
					continue
				else:
					if flag == 0:
						output = self.module[i]+" total hits = "+str(self.hits[i])
						print output
						flag += 1
					if len(self.sigRes[i][j]) == 0:
						print "The ",j+1,"th has no significant dis"
					else:
						output = "The "+str(j+1)+"th significant dis: "+str(self.sigRes[i][j])
						print output
						
	def displayPWMresult(self,thresHits,motif_obj):
		for i in range(len(motif_obj.moduleName)):
			flag = 0
			for j in range(len(self.sigRes[i])):
				if self.sigRes[i] == "null" or self.hits[i] < thresHits:
					continue
				else:
					if flag == 0:
						output = motif_obj.moduleName[i]+" total hits = "+str(self.hits[i])
						print output
						flag += 1
					if len(self.sigRes[i][j]) == 0:
						print "The ",j+1,"th has no significant dis"
					else:
						output = "The "+str(j+1)+"th significant dis: "+str(self.sigRes[i][j])
						print output
						
	def displayNewModule(self,module_array):
		output_module=[]
		suffix = ".modres"
		for module in module_array:
			f = open('pattern_file','w')
			f.write(module)
			f.close()
			result_file = module + suffix
			command = PatScan+"pattern_file < " + self.Promoter
			lines = os.popen(command).readlines()
			self.dis_result.append(lines)
			total_hits = len(lines)/2
			if total_hits < thresHits:
				continue
			output_module.append(module)
			f = open(result_file,'w')
			f.write(module)
			f.write("\n")
			for line in lines:
				f.write(line)
				f.write("\n")
			print module," total hits: ",total_hits
		if len(module_array) <= 1:
			return [0,output_module]
		else:
			return [1,output_module]		
		