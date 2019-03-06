class Prep():
    def __init__(self):
        self.ndict = {}
        self.edict = {}
        self.nnum = 0
        self.enum = 0

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

    def nswap_in_eles(self, nnum_list, nnum_tomerge):
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

    def nclear_unused(self, nnum_list, nnum_tomerge):
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
        return self.supnum

    def kpoint(self, x, y):
        self.knum += 1
        self.kdict[self.knum] = [x, y]
        return self.kdict[self.knum]

class Solv():
    def __init__(self, meshclass):
        self.prepclass = meshclass.PREP
        self.meshclass = meshclass
        nodesnum = len(self.prepclass.ndict)
        self.clist = [0 for i in range(nodesnum * 2)]
        self.gklist = [[0 for i in range(nodesnum * 2)] for j in range(nodesnum * 2)]

    @staticmethod
    def stiff_matrix(E, v, h):
        """Element stiffness matrix"""
        p = (E * h) / (12 * (1 - (v ** 2)))
        q = (E * h) / (1 - v)
        k11 = p * 2 * (3 - v)
        k15 = p * (-3 + v)
        k12 = q / 8
        k16 = -k12
        k17 = p * 2 * v
        k18 = p * 1.5 * (1 - (3 * v))
        k13 = p * (-3 - v)
        k14 = -k18
        klist = [[k11, k12, k13, k14, k15, k16, k17, k18],
                 [k12, k11, k18, k17, k16, k15, k14, k13],
                 [k13, k18, k11, k16, k17, k14, k15, k12],
                 [k14, k17, k16, k11, k18, k13, k12, k15],
                 [k15, k16, k17, k18, k11, k12, k13, k14],
                 [k16, k15, k14, k13, k12, k11, k18, k17],
                 [k17, k14, k15, k12, k13, k18, k11, k16],
                 [k18, k13, k12, k15, k14, k17, k16, k11]]
        return klist
    
    def build(self):
        for meshnum in self.meshclass.meshdict:
            eles = self.meshclass.meshdict[meshnum][0]
            mat_params = self.meshclass.meshdict[meshnum][1]
            thickness = self.meshclass.meshdict[meshnum][2]
            size = self.meshclass.meshdict[meshnum][3]
            
            klist = self.stiff_matrix(mat_params[0], mat_params[1], thickness)
            
            for ele in eles:
                pass
    
class Mesh():
    def __init__(self, geomclass):
        self.geomclass = geomclass
        self.PREP = Prep()
        self.meshdict = {}
        self.meshnum = 0
        self.next_nnum = 1
        
    def generate(self, rectnum, size=1):
        self.meshnum += 1
        nodes = []
        eles = []
        
        rectgeom = self.geomclass.supdict[rectnum]
        xlist, ylist = [rectgeom[0][0], rectgeom[1][0]], [rectgeom[0][1], rectgeom[1][1]]
        width = abs(rectgeom[1][0] - rectgeom[0][0])
        height = abs(rectgeom[1][1] - rectgeom[0][1])

        for i in range(min(xlist), max(xlist) + 1, size):
            for j in range(min(ylist), max(ylist) + 1, size):
                nodes.append(self.PREP.node(i, j))

        for i in range(width):
            for j in range(self.next_nnum, self.next_nnum + height):
                eles.append(self.PREP.ele(j + 1 + ((height + 1) * i), 
                                          j + 1 + ((height + 1) * (i + 1)),
                                          j + ((height + 1) * (i + 1)), 
                                          j + ((height + 1) * i)))

        self.meshdict[self.meshnum] = [eles, 0, 0, size]
        self.next_nnum += 1 + nodes[-1] - nodes[0]
        return self.meshnum

    def nmerge(self, x, y):
        nnum_list = self.PREP.ncheck(x, y)
        if nnum_list:
            nnum_tomerge = min(nnum_list)
            self.PREP.nswap_in_eles(nnum_list, nnum_tomerge)
            self.PREP.nclear_unused(nnum_list, nnum_tomerge)

    def assignprop(self, meshnum, mat_params, thickness):
        self.meshdict[meshnum][1] = mat_params
        self.meshdict[meshnum][2] = thickness