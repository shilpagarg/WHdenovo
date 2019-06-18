import sys
mainset = set(open(sys.argv[1],'r').readlines())
exc = set(open(sys.argv[2], 'r').readlines())
diff = mainset.difference(exc)
o = open(sys.argv[3], 'w')
for i in diff:
    o.write(i)
o.close()
