import re

def readkey(filename):

    prm_sigma = {}
    prm_ep = {}
    prm_ch = {}

    keyfile = open(filename,'r')

    while True:
        templine = keyfile.readline()
        if templine == '': break
        line = templine.strip()
        linedata = re.split('\s+', line)
        if linedata[0] == 'parameters': prmfilename = linedata[1]
        elif linedata[0] == 'vdw':
            prm_sigma[int(linedata[1])] = float(linedata[2])
            prm_ep[int(linedata[1])] = float(linedata[3])
        elif linedata[0] == 'charge': prm_ch[int(linedata[1])] = float(linedata[2])

    keyfile.close()

    prmfile = open(prmfilename,'r')

    while True:
        templine = prmfile.readline()
        if templine == '': break
        line = templine.strip()
        linedata = re.split('\s+', line)
        if linedata[0] == 'electric': coul_conv = float(linedata[1])
        if linedata[0] == 'vdw':
            prm_sigma[int(linedata[1])] = float(linedata[2])
            prm_ep[int(linedata[1])] = float(linedata[3])
        elif linedata[0] == 'charge': prm_ch[int(linedata[1])] = float(linedata[2])

    prmfile.close()

    output = [prm_sigma,prm_ep,prm_ch,coul_conv]

    return output
