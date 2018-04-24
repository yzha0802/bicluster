def sumNums(numbers):
    total=0
    for num in numbers:
            total+=num
    return total

def averageNums(numbers):
        sumOfNums =sumNums(numbers)
        average=float(sumOfNums) / len(numbers)
        return average
def varianceNums(numbers):
        variance =[0]*len(numbers)
        average = averageNums(numbers)
        for num in numbers:
                variance[numbers.index(num)]=(num-average)**2
        return averageNums(variance)
def stdDevNums(numbers):
        variance= varianceNums(numbers)
        try:
            return variance ** .5
        except (TypeError):
                print("wrong data")
