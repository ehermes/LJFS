import re

def readarc(filename):

    arcfile = open(filename,'r')

    n = 0
    atomname = []
    x = []
    y = []
    z = []
    atomtype = []
    bonds = []

    while True:
        atomnamen = []
        xn = []
        yn = []
        zn = []
        atomtypen = []
        bondsn = []
        line = arcfile.readline().strip()
        if line == '': break
        data = re.split('\s+', line)
        count = int(data[0])
        for i in xrange(1, count+1):
            line = arcfile.readline().strip()
            data = re.split('\s+', line)
            number = int(data[0])
            assert number == i
            atomnamen.append(data[1])
            xn.append(float(data[2]))
            yn.append(float(data[3]))
            zn.append(float(data[4]))
            atomtypen.append(int(data[5]))
            bondsn.append([int(k) for k in data[6:]])
        atomname.append(atomnamen)
        x.append(xn)
        y.append(yn)
        z.append(zn)
        atomtype.append(atomtypen)
        bonds.append(bondsn)
        n += 1

    output = [n, atomname, x, y, z, atomtype, bonds]

    return output
