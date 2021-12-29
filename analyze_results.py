import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

results = []  # in the form aw ah bw bh iterations num_threads mean

mypath = "results/"
for r, d, f in os.walk(mypath):
    for file in f:
        if file.endswith(".txt"):
            txt = os.path.join(r, file).split('.')[0].split('_')
            aw = int(txt[1].split('x')[0])
            ah = int(txt[1].split('x')[1])
            bw = int(txt[2].split('x')[0])
            bh = int(txt[2].split('x')[1].split('/')[0])
            algorithm = 0
            if txt[4] != "time":
                iterations = int(txt[4].split('t')[1])
                if len(txt) > 5:
                    num_threads = int(txt[5].split('h')[1])
                else:
                    num_threads = 1
            else: # regular
                algorithm = 1
                iterations = int(txt[5].split('t')[1])
                if len(txt) > 6:
                    num_threads = int(txt[6].split('h')[1])
                else:
                    num_threads = 1

            curr_times = []
            time_file = open(os.path.join(r, file), "r")
            for line in time_file:
                curr_times.append(list(map(float, line.split())))
            curr_times = np.matrix(curr_times)
            mean = curr_times.mean()

            results.append([ah, aw, bw, bh, iterations, num_threads, mean, algorithm])

# print(results)
threads = []
for i in results:
    if i[5] not in threads:
        threads.append(i[5])
threads.sort()

mat_size = []
for i in results:
    if i[0] not in mat_size:
        mat_size.append(i[0])
mat_size.sort()

for i in mat_size:
    line = np.zeros(len(threads))
    for j in results:
        if j[0] == i:
            if j[7] == 0:
              line[j[5]-1] = j[6]
    plt.plot(threads, line)

    # line = np.zeros(len(threads))
    # for j in results:
    #     if j[0] == i:
    #         if j[7] == 1:
    #           line[j[5]-1] = j[6]
    # plt.plot(threads, line)

    plt.xlabel("Number of threads")
    plt.ylabel("Average time [ns]")
    plt.title("Matrix size: " + str(i))
    plt.savefig("results/mat_size_" + str(i) + ".png")
    # plt.figure()
