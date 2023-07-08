import scipy.io as sio
import sys

mat = sio.loadmat(sys.argv[1])
data = mat["G"][0][0]
print("ID\tn\tp\tq")
d = []
for i in range(len(data[0])):
    #print("%s\t%s\t%s\t%s" % (data[0][i][0][0], data[72][i][0], data[41][i][0], data[66][i][0]) )
    d.append([data[0][i][0][0], str(data[72][i][0]), str(data[41][i][0]), str(data[66][i][0])])

d = sorted(d, key=lambda x: float(x[3]))
for l in d:
    print("\t".join(l))
