import sys
a = set(open(sys.argv[1], 'r').readlines())
sys.stderr.write('uniq items in first file: %d\n'%len(a))
b = set(open(sys.argv[2], 'r').readlines())
sys.stderr.write('uniq items in first file: %d\n'%len(b))
c = a.intersection(b)
sys.stderr.write('size of the intersection: %d\n'%len(c))
for i in c:
    sys.stdout.write(i)

