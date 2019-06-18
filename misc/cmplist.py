import subprocess
total = 0
tw = 0
for i in range(569):

    diff1 = subprocess.Popen('diff -y block_%d ~/WHdenovo/trioasm/whatshap_trioasm/refBasedBC.ordered.nodes| grep -v \'>\''%i, shell = True, stderr = subprocess.PIPE, 
                            stdout = subprocess.PIPE)
    com1 = diff1.communicate()
    res1 = com1[0].decode().split('\n')

    diff2 = subprocess.Popen('diff -y block_%d_rev ~/WHdenovo/trioasm/whatshap_trioasm/refBasedBC.ordered.nodes| grep -v \'>\''%i, shell = True, stderr = subprocess.PIPE, 
                            stdout = subprocess.PIPE)
    com2 = diff2.communicate()
    res2 = com2[0].decode().split('\n')

    w1 = 0
    for line in res1:
        if '|' in line or '<' in line:
            w1 += 1
    if w1 == 0:
        total += len(res1)
        continue

    w2 = 0
    for line in res2:
        if '|' in line or '<' in line:
            w2 += 1
    if w2 == 0:
        total += len(res2)
        continue    

    print('block_%d'%i)
    if w1 <= w2:
        for b in res1:
            print(b)
    else:
        for b in res2:
            print(b)

    total += min(len(res1), len(res2))
    tw += min(w1, w2)

print('total nodes/bubbles from algorithm', total)
print('wrong ones', tw)
