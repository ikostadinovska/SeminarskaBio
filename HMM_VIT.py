__author__ = 'Ivana'

import math

"""  function for creating the hidden markov model"""
def makeHMM(s, h):
    t = {}
    e = {}
    prev_hs = 0
    print zip(h,s)
    for hs, ss in zip(h, s):

        tmpd = e.setdefault(hs, {})
        tmpd[ss] = tmpd.get(ss, 0) + 1
        tmpd = t.setdefault(prev_hs, {})
        tmpd[hs] = tmpd.get(hs, 0) + 1

        prev_hs = hs
    nt = {}
    for h1, h2_f in t.iteritems():
        n = float(sum(h2_f.values()))
        tmpd = {}
        for h2, f in h2_f.iteritems():
            rf = f / n
            tmpd[h2] = rf
        nt[h1] = tmpd
    ne = {}
    for h1, s1_f in e.iteritems():
        n = float(sum(s1_f.values()))
        tmpd = {}
        for s1, f in s1_f.iteritems():
            rf = f / n
            tmpd[s1] = rf
        ne[h1] = tmpd
    hmm = (nt, ne)
    return hmm


"""  function for implementing the Viterbi algorithm"""
def Viterbialgorithm(s, hmm):
    min_val = 0.000001
    t, e = hmm

    lh = set()
    for h, tmpd in e.iteritems():
        lh.add(h)
    lh = [0] + list(lh)

    tableV = [{} for i in range(len(s)+1)]
    ptr = [{} for i in range(len(s)+1)]

    for k in lh:
        tableV[0][k] = math.log(min_val)
    tableV[0][0] = math.log(1.0)

    for i in range(1, len(s)+1):
        for l in lh:
            vals = [(tableV[i-1][k]+math.log(t[k].get(l, min_val)), k) for k in lh]
            max_val, max_k = max(vals)
            tableV[i][l] = math.log(e.get(l, {}).get(s[i-1], min_val)) + max_val
            ptr[i][l] = max_k

    pi = []
    pi_L = max([(tableV[-1][k], k) for k in lh])[1]
    pi.append(pi_L)

    for p in ptr[-1:1:-1]:
        pi.append(p[pi[-1]])
    pi.reverse()
    return "".join(pi)



###################   making  test data, hidden and observed sequence for k=1
f = open("sequence_sacCer_training.fasta")
fc = f.read()
seq = fc.split("\n")
seqObserved = list(seq[1]) # Observed sequence for k=1

ft = open("sequence_sacCer_test.fasta")
fct = ft.read()
seq2 = fct.split("\n")
testData = list(seq2[1])  # test data a for k=1

# Hidden sequence
seqHidden = []
listTraining = open("annotation_sacCer_training.txt", "rt").readlines()

for i in range(1,len(listTraining)):
    seqHidden.append(listTraining[i][6])


###############################  KREIRANJE NA K-MERS = 2  #################################

# create seqObserved k=2
seqPOMOS = seq[1]
seqOb_2 = []
for i in range(len(seqPOMOS)-1):
    seqOb_2.append(seqPOMOS[i]+seqPOMOS[i+1])
seqOb_2.append(seqPOMOS[len(seqPOMOS)-1]+"A")

#create test data for k=2
seqPOM = seq2[1]
seqTest_2 = []
for i in range(len(seqPOM)-1):
    seqTest_2.append(seqPOM[i]+seqPOM[i+1])
seqTest_2.append(seqPOM[len(seqPOM)-1]+"A")


print seqOb_2[:12]
print seqTest_2[:12]
#call function for creating HMM
trainLength = 96598 #65%
HMM = makeHMM(seqOb_2[:trainLength/2], seqHidden[:trainLength/2])
resultH = Viterbialgorithm(seqTest_2, HMM)



results = open("results_VIT_4_k_2.txt", "w")
for i in range(0, len(resultH)):
    results.write(resultH[i] + '\n')
results.close()

print "Hidden Markov Model: ", HMM
#print "hidden state sequence:", resultH



#
########### VIZUELIZACILA
#
#
#
#import matplotlib.pyplot as plt
#x_oska = testData[74900:75000] #nucleotidite  143500
#y_oska = resultH[74900:75000] #hidden states
#
#x = [i for i in range(len(x_oska)-1)]
#y = [1 if(y_oska[i] == "C") else 0 for i in range(len(y_oska)-1) ]
#labels = [x_oska[i] for i in range(len(x_oska)-1)]
#
#plt.plot(x, y, 'ro')
## You can specify a rotation for the tick labels in degrees or with keywords.
#plt.xticks(x, labels, rotation='horizontal')
## Pad margins so that markers don't get clipped by the axes
#plt.margins(0.2)
## Tweak spacing to prevent clipping of tick-labels
#plt.subplots_adjust(bottom=0.15)
#plt.show()