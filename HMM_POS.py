__author__ = 'Ivana'

import math

"""  function for creating the hidden markov model"""
def makeHMM(s, h):
    t = {}
    e = {}
    prev_hs = 0
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



def logmv(a):
    min_val = 0.0000000001
    return math.log(max(a, min_val))

def log_sum(v1, v2):
    nv1 = max(v1, v2)
    nv2 = min(v1, v2)
    return nv1 + math.log(1.0+math.exp(nv2-nv1))

def sum_log(vals):
    s = vals[0]
    for v in vals[1:]:
        s = log_sum(s, v)
    return s

"""  function for creating the forward and backward algorithms"""
def for_back(s, hmm):
    t, e = hmm
    zh = e.keys()
    zh = [0] + list(zh)
    f = [{} for i in range(len(s)+1)]
    b = [{} for i in range(len(s)+1)]
    for k in zh:
        f[0][k] = logmv(0.0)
        b[len(s)][k] = math.log(1.0)
    f[0][0] = math.log(1.0)

    for i in range(1, len(s)+1):
        for l in zh:
            sum_val = sum_log([f[i-1][k] + logmv(t[k].get(l, 0.0)) for k in zh])
            f[i][l] = logmv(e.get(l, {}).get(s[i-1], 0.0)) + sum_val
    ps_for = sum_log([f[len(s)][k] for k in zh])

    for i in range(len(s)-1,0,-1):
        for k in zh:
            sum_val = sum_log([logmv(t[k].get(l, 0.0)) + logmv(e.get(l, {}).get(s[i], 0.0)) + b[i+1][l] for l in zh])
            b[i][k] = sum_val
    ps_back = sum_log([logmv(t[0].get(l, 0.0)) + logmv(e.get(l, {}).get(s[0], 0.0)) + b[1][l] for l in zh])

    return f, ps_for, b, ps_back

"""  function for creating the posterior decoding algorithm"""
def posterior_decoding(s, hmm):
    f, fps, b, bps = for_back(s, hmm)
    if not (abs(fps - bps) < 0.01):
        print fps, bps, fps-bps
        print "ERROR"
        return ""

    hn = ""
    for i in range(1, len(s)+1):
        max_p, max_k = max([(fp+b[i][k], k) for k, fp in f[i].iteritems()])
        hn = hn + max_k
    return hn




#making  test data, hidden and observed sequence for k=1
f = open("sequence_sacCer_training.fasta")
fc = f.read()
seq = fc.split("\n")
seqObserved = list(seq[1])  # Observed sequence for k=1

ft = open("sequence_sacCer_test.fasta")
fct = ft.read()
seq2 = fct.split("\n")
testData = list(seq2[1]) #test data a for k=1

# Hidden sequence
seqHidden = []
listTraining = open("annotation_sacCer_training.txt", "rt").readlines()

for i in range(1,len(listTraining)):
    seqHidden.append(listTraining[i][6])



#call function for creating HMM
trainLength = len(seqObserved) #65%
HMM = makeHMM(seqObserved[:trainLength/2], seqHidden[:trainLength/2])
resultH = posterior_decoding(testData, HMM)

results = open("results_POS_5.txt", "w")
for i in range(0, len(resultH)):
    results.write(resultH[i] + '\n')
results.close()

print "Hidden Markov Model: ", HMM
print "hidden state sequence:", resultH
