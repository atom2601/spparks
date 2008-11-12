# Pizza.py script to convert a SPPARKS dump file to a SPPARKS input file
# assumes dump file is in LAMMPS dump format
# input file is used by lattice apps with "input" keyword in app_style

# Syntax: pizza.py -f dump2init.py dumpfile N inputfile
# N = timestamp of snapshot to be converted

if len(argv) != 4:
  raise StandardError, "Syntax: pizza.py -f dump2init.py dumpfile N inputfile"

dumpfile = argv[1]
ntime = int(argv[2])
inputfile = argv[3]

# edit number of columns and column assignments if needed

d = dump(dumpfile)
d.map(1,"id",2,"lattice")
id,lattice = d.vecs(ntime,"id","lattice")

fp = open(inputfile,"w")
print >>fp,"SPPARKS input file from dump file, time =",ntime
print >>fp,len(id),1
print >>fp
for i in xrange(len(id)):
  print >>fp,int(id[i]),int(lattice[i])
fp.close()
