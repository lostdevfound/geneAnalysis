import numpy as np
import sys
from matplotlib import pyplot as plt
import seaborn as sb

if __name__ == "__main__":
    k = sys.argv[1]
    fname1 = sys.argv[2]
    data1 = np.loadtxt(fname1, dtype='int', usecols=[1])

    plt.plot(data1, label=fname1)
    plt.title("Coverage histogram, k={}".format(k))
    plt.axis([1,500, 1, 500])
    plt.ylabel('number of k-mers')
    plt.xlabel('coverage')
    plt.xticks(np.arange(0,500,10))     # set tick freq
    plt.yticks(np.arange(0,1000,100))     # set tick freq


    if len(sys.argv) > 3:
        fname2 = sys.argv[3]
        data2 = np.loadtxt(fname2, dtype='int', usecols=[1])
        plt.plot(data2, 'g', label=fname2)
    plt.show()
