import module_detection3
import sys
from time import time
import os

start = time()
test = module_detection3.motifs(sys.argv[1])
test.readFile()

module = test.motif
flag = 1
while flag != 0:
	#print "Assembly modules"
	new = test.assemblyModules(module)
	a = module_detection3.modules(new)
	#print "Start searching Modules..."
	a.searchModule(sys.argv[2])
	#print "Testing Statistics..."
	module = a.testZScore()
	#print "Displaying modules"
	[flag,module] = a.displayNewModule(module)

end = time()
print "total time is: "
print end-start

if os.path.exists("result") == False:
	os.mkdir("result")
os.popen("mv *.modres result")

