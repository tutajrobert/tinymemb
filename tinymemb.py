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
            for i in range(4):
                if self.edict[enum][i] in nnum_list:
                    self.edict[enum][i] = nnum_tomerge

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
        self.nod2dofmap = {}
        self.dofnum = 0
        self.nodesnum = len(self.prepclass.ndict)
        self.clist = [0 for i in range(self.nodesnum * 2)]
        self.gklist = [[0 for i in range(self.nodesnum * 2)] for j in range(self.nodesnum * 2)]
        self.dofmapping()

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

    def dofmapping(self):
        for nnum in self.prepclass.ndict:
            self.nod2dofmap[nnum] = self.dofnum
            self.dofnum += 2

    def build(self):
        for meshnum in self.meshclass.meshdict:
            eles = self.meshclass.meshdict[meshnum][0]
            mat_params = self.meshclass.meshdict[meshnum][1]
            thickness = self.meshclass.meshdict[meshnum][2]
            size = self.meshclass.meshdict[meshnum][3]
            klist = self.stiff_matrix(mat_params[0], mat_params[1], thickness)
            for enum in eles:
                nodes = self.prepclass.edict[enum]
                dofs = []
                for nnum in nodes:
                    dofs.append(self.nod2dofmap[nnum])
                    dofs.append(self.nod2dofmap[nnum] + 1)
                for i in range(8):
                    for j in range(8):
                        self.gklist[dofs[i]][dofs[j]] += klist[i][j]

    def support(self, x, y):
        nnum = self.prepclass.ncheck(x, y)[0]
        dofs = [self.nod2dofmap[nnum], self.nod2dofmap[nnum] + 1]
        for dof in dofs:
            self.gklist[dof] = [0 for i in range(self.nodesnum * 2)]
            for row in self.gklist:
                row[dof] = 0
            self.gklist[dof][dof] = 1

    def force(self, x, y, force_value):
        nnum = self.prepclass.ncheck(x, y)[0]
        dof = self.nod2dofmap[nnum]
        self.clist[dof] = force_value
    
    def connector(self, node1, node2):
        dofs = []
        dof1 = self.nod2dofmap[node1]
        dof2 = self.nod2dofmap[node2]
        dofs.append(dof1)
        dofs.append(dof1 + 1)
        dofs.append(dof2)
        dofs.append(dof2 + 1)
        k = 2e10
        klist = [[k, 0, -k, 0],
                 [0, k, 0, -k],
                 [-k, 0, k, 0],
                 [0, -k, 0, k]]
        
        for i in range(4):
            for j in range(4):
                self.gklist[dofs[i]][dofs[j]] += klist[i][j]

    @staticmethod
    def gausselim(m):
        #eliminate columns
        for col in range(len(m[0])):
            for row in range(col+1, len(m)):
                r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
                m[row] = [sum(pair) for pair in zip(m[row], r)]
        #now backsolve by substitution
        ans = []
        m.reverse() #makes it easier to backsolve
        for sol in range(len(m)):
            if sol == 0:
                ans.append(m[sol][-1] / m[sol][-2])
            else:
                inner = 0
                #substitute in all known coefficients
                for x in range(sol):
                    inner += (ans[x]*m[sol][-2-x])
                #the equation is now reduced to ax + b = c form
                #solve with (c - b) / a
                ans.append((m[sol][-1]-inner)/m[sol][-sol-2])
        ans.reverse()
        return ans

    def solve(self):
        for i in range(len(self.clist)):
            self.gklist[i].append(self.clist[i])
        self.dlist = self.gausselim(self.gklist)
        return self.dlist

    def backsolve4force(self, x, y, dlist):
        nnum = self.prepclass.ncheck(x, y)[0]
        print(nnum)
        dof = self.nod2dofmap[nnum]
        force = 0
        for i in range(0, self.nodesnum * 2):
            if self.gklist[dof][i] * dlist[i]:
                force += self.gklist[dof][i] * dlist[i]
        return force

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
        width = int(abs(rectgeom[1][0] - rectgeom[0][0]) / size)
        height = int(abs(rectgeom[1][1] - rectgeom[0][1]) / size)

        for i in range(min(xlist), max(xlist) + 1, size):
            for j in range(min(ylist), max(ylist) + 1, size):
                nodes.append(self.PREP.node(i, j))

        for i in range(width):
            for j in range(self.next_nnum, self.next_nnum + height):
                eles.append(self.PREP.ele(j + 1 + ((height + 1) * i), 
                                          j + 1 + ((height + 1) * (i + 1)),
                                          j + ((height + 1) * (i + 1)), 
                                          j + ((height + 1) * i)))
        self.meshdict[self.meshnum] = [eles, 0, 0, size, nodes]
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
        
    def selectnode(self, coords, meshnum):
        x, y = coords[0], coords[1]
        nnum_list = self.PREP.ncheck(x, y)
        for nnum in nnum_list:
            if nnum in self.meshdict[meshnum][4]:
                return nnum