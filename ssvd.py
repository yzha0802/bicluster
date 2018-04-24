
myFile=open("example.fastq","r")
adapterSequence='GCCAAT'
totalLines=0
countOfAdapter=0
for line in myFile:
        if line[0]=='N':
                if adapterSequence in line:
                        countOfAdapter+=1
                totalLines+=1
print("Total Lines:%.0f" % totalLines)
print("Count of adapter:%.0f") % countOfAdapter
percentage=(float(countOfAdapter)/totalLines)*100
print("Percentage of reads containg the adapter:%2f" % percentage)
