class Prep():
    def __init__(self):
        self.ndict = {}
        self.edict = {}
        self.nnum = 0
        self.enum = 0
        self.meshnum = 0
        self.meshdict = {}
        self.next_nnum = 1

    def node(self, x, y):
        self.nnum += 1
        self.ndict[self.nnum] = [x, y]
        return self.nnum

    def ndel(self, nnum):
        del self.ndict[nnum]

    def ncheck(self, x, y):
        nnum_list = []
        for nnum in self.ndict:
            if self.ndict[nnum] == [x, y]:
                nnum_list.append(nnum)
        return nnum_list

    def ele(self, n1, n2, n3, n4):
        self.enum += 1
        self.edict[self.enum] = [n1, n2, n3, n4]
        return self.enum


    def mesh(self, supele, size, property):
        self.meshnum += 1
        nodes = []
        eles = []

        start, end = supele[0], supele[1]
        xlist, ylist = [start[0], end[0]], [start[1], end[1]]
        xmin, xmax = min(xlist), max(xlist)
        ymin, ymax = min(ylist), max(ylist)
        xlen = abs(end[0] - start[0])
        ylen = abs(end[1] - start[1])

        for i in range(xmin, xmax + 1, size):
            for j in range(ymin, ymax + 1, size):
                nodes.append(self.node(i, j))

        for i in range(xlen):
            for j in range(self.next_nnum, self.next_nnum + ylen):
                eles.append(self.ele(j + 1 + ((ylen + 1) * i), j + 1 + ((ylen + 1) * (i + 1)),
                                     j + ((ylen + 1) * (i + 1)), j + ((ylen + 1) * i)))

        self.meshdict[self.meshnum] = [eles]
        self.next_nnum += 1 + nodes[-1] - nodes[0]
        return self.meshnum

    def nmerge(self, x, y):
        nnum_list = self.ncheck(x, y)
        if len(nnum_list) > 1:
            nnum_tomerge = min(nnum_list)
        for enum in self.edict:
            for nnum in nnum_list:
                if self.edict[enum][0] == nnum:
                    self.edict[enum][0] = nnum_tomerge
                elif self.edict[enum][1] == nnum:
                    self.edict[enum][1] = nnum_tomerge
                elif self.edict[enum][2] == nnum:
                    self.edict[enum][2] = nnum_tomerge
                elif self.edict[enum][3] == nnum:
                    self.edict[enum][3] = nnum_tomerge
        for nnum in nnum_list:
            if nnum != nnum_tomerge:
                self.ndel(nnum)

class Geom():
    def __init__(self):
        self.kdict = {}
        self.supdict = {}
        self.supnum = 0
        self.knum = 0

    def supele(self, kpoint_start, kpoint_end):
        self.supnum += 1
        self.supdict[self.supnum] = [kpoint_start, kpoint_end]
        return self.supdict[self.supnum]

    def kpoint(self, x, y):
        self.knum += 1
        self.kdict[self.knum] = [x, y]
        return self.kdict[self.knum]

class Solv():
    def __init__(self, prepclass):
        self.prepclass = prepclass
        nodesnum = len(self.prepclass.ndict)
        self.clist = [0 for i in range(nodesnum * 2)]
        self.gklist = [[0 for i in range(nodesnum * 2)] for j in range(nodesnum * 2)]

    @staticmethod
    def stiff_matrix():
        """Element stiffness matrix"""
        v = .5
        p = 1 / (12 * (1 - (v ** 2)))
        q = 1 / (1 - v)
        k11 = p * 2 * (3 - v)
        k15 = p * (-3 + v)
        k12 = q / 8
        k16 = -k12
        k17 = p * 2 * v
        k18 = p * 1.5 * (1 - (3 * v))
        k13 = p * (-3 - v)
        k14 = -k18

        #Element stiffnes matrix
        klist = [[k11, k12, k13, k14, k15, k16, k17, k18],
                 [k12, k11, k18, k17, k16, k15, k14, k13],
                 [k13, k18, k11, k16, k17, k14, k15, k12],
                 [k14, k17, k16, k11, k18, k13, k12, k15],
                 [k15, k16, k17, k18, k11, k12, k13, k14],
                 [k16, k15, k14, k13, k12, k11, k18, k17],
                 [k17, k14, k15, k12, k13, k18, k11, k16],
                 [k18, k13, k12, k15, k14, k17, k16, k11]]
        return klist