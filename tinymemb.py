"""Tiny finite element method application for plane stress problems and bolt's forces"""
class Prep():
    """Preprocessor class. Not used by the user"""
    def __init__(self):
        self.ndict = {}
        self.edict = {}
        self.nnum = 0
        self.enum = 0

    def node(self, x, y):
        """Creates node and returns its number"""
        self.nnum += 1
        self.ndict[self.nnum] = [x, y]
        return self.nnum

    def deletenode(self, nodenum):
        """Deletes node entry from nodes dictionary"""
        del self.ndict[nodenum]

    def nodes_bycoords(self, x, y):
        """Returns node numbers of given x, y position in list"""
        nodenum_list = []
        for nodenum in self.ndict:
            if self.ndict[nodenum] == [x, y]:
                nodenum_list.append(nodenum)
        return nodenum_list

    def element(self, n1, n2, n3, n4):
        """Creates element and returns its number"""
        self.enum += 1
        self.edict[self.enum] = [n1, n2, n3, n4]
        return self.enum

    def swap_nodenum(self, nodenum_toleft, nodenum_tomerge):
        """Swaps node numbers to merge two nodes"""
        for enum in self.edict:
            for i in range(4):
                if self.edict[enum][i] == nodenum_tomerge:
                    self.edict[enum][i] = nodenum_toleft

class Geom():
    """Geometry class contains points and rectangles"""
    def __init__(self):
        self.rectsdict, self.pointsdict = {}, {}
        self.rectnum, self.pointnum = 0, 0

    def rectangle(self, startpoint, endpoint):
        """Creates rectangle from startpoint to endpoint. Returns rectangle number"""
        self.rectnum += 1
        self.rectsdict[self.rectnum] = [startpoint, endpoint]
        return self.rectnum

    def point(self, x, y):
        """Creates point of x, y coordinates. Returns coordinates list [x, y]"""
        self.pointnum += 1
        self.pointsdict[self.pointnum] = [x, y]
        return self.pointsdict[self.pointnum]

class Solv():
    """Solver class contains matrices and functions needed for calculation"""
    def __init__(self, meshclass):
        self.prepclass = meshclass.prepclass
        self.meshclass = meshclass
        self.nodesnum = len(self.prepclass.ndict) #number of nodes
        self.fvector = [0 for i in range(self.nodesnum * 2)] #external forces vector
        self.kmatrix = [[0 for i in range(self.nodesnum * 2)] for j in range(self.nodesnum * 2)]
        self.kmatrix_4backsolve = [[0 for i in range(self.nodesnum * 2)] for j in range(self.nodesnum * 2)]
        self.nodenum_todof_mapping()

    @staticmethod
    def stiffmatrix_ofsquare(E, v, h):
        """Returns 4 noded membrane element stiffness matrix"""
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

    def nodenum_todof_mapping(self):
        """Creates node number to degree of freedom relation"""
        self.nodenum_todof = {} #dictionary for nodes to dofs mapping
        self.dofnum = 0
        for nodenum in self.prepclass.ndict:
            self.nodenum_todof[nodenum] = self.dofnum
            self.dofnum += 2

    def build(self):
        """Global stiffness matrix aggregation"""
        for meshnum in self.meshclass.meshdict:
            eles = self.meshclass.meshdict[meshnum][0]
            mat_params = self.meshclass.meshdict[meshnum][1] #[young modulus, poisson ratio]
            thickness = self.meshclass.meshdict[meshnum][2]
            elematrix = self.stiffmatrix_ofsquare(mat_params[0], mat_params[1], thickness)
            for enum in eles:
                nodes = self.prepclass.edict[enum]
                dof_list = []
                for nnum in nodes:
                    dof_list.append(self.nodenum_todof[nnum])
                    dof_list.append(self.nodenum_todof[nnum] + 1)
                for i in range(8):
                    for j in range(8):
                        self.kmatrix[dof_list[i]][dof_list[j]] += elematrix[i][j]
                        self.kmatrix_4backsolve[dof_list[i]][dof_list[j]] += elematrix[i][j]
    def support(self, nodenum, x=True, y=True):
        """Supports creation procedure"""
        dof_list = [self.nodenum_todof[nodenum], self.nodenum_todof[nodenum] + 1]
        if x is False:
            dof_list.pop(0)
        if y is False:
            dof_list.pop(1)
        for dof in dof_list:
            self.kmatrix[dof] = [0 for i in range(self.nodesnum * 2)]
            for row in self.kmatrix:
                row[dof] = 0
            self.kmatrix[dof][dof] = 1

    def force(self, nodenum, xforce=0, yforce=0):
        """Force application procedure"""
        dof = self.nodenum_todof[nodenum]
        self.fvector[dof] = xforce
        self.fvector[dof + 1] = yforce

    def connector(self, node1, node2, k=2e10):
        """Inplane onnector creation for nodes pair"""
        dofs = []
        dof1 = self.nodenum_todof[node1]
        dof2 = self.nodenum_todof[node2]
        dofs.append(dof1)
        dofs.append(dof1 + 1)
        dofs.append(dof2)
        dofs.append(dof2 + 1)
        connmatrix = [[k, 0, -k, 0],
                      [0, k, 0, -k],
                      [-k, 0, k, 0],
                      [0, -k, 0, k]]
        for i in range(4):
            for j in range(4):
                self.kmatrix[dofs[i]][dofs[j]] += connmatrix[i][j]

    @staticmethod
    def gausselim(m):
        """Gauss elimination procedure"""
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
                #the equation is now reduced to ax + b = c
                #solve with (c - b) / a
                ans.append((m[sol][-1]-inner)/m[sol][-sol-2])
        ans.reverse()
        return ans

    def solve(self):
        """Solve matrix equation kmatrix * uvector = fvector"""
        for i in range(len(self.fvector)):
            self.kmatrix[i].append(self.fvector[i])
        self.uvector = self.gausselim(self.kmatrix)
        return self.uvector

    def backsolve_4force(self, nodenum):
        """Back calculation for force"""
        dof = self.nodenum_todof[nodenum]
        force = 0
        for i in range(0, self.nodesnum * 2):
            if self.kmatrix_4backsolve[dof][i] * self.uvector[i]:
                force += self.kmatrix_4backsolve[dof][i] * self.uvector[i]
        return force

class Mesh():
    """Mesh class contains meshes and initializes preprocessor class"""
    def __init__(self, geomclass):
        self.geomclass = geomclass #access to rectangles needed
        self.prepclass = Prep() #initialization of preprocessor class
        self.meshdict = {}
        self.meshnum = 0
        self.nodenum_4nextmesh = 1 #number of next node, which is a starting number for every mesh

    def generate(self, rectnum, size=1):
        """Mesh generation procedure. Creates nodes and elements in preprocessor class.
           Returns generated mesh number"""
        self.meshnum += 1
        nodes = []
        eles = []
        rectgeom = self.geomclass.rectsdict[rectnum]
        xlist, ylist = [rectgeom[0][0], rectgeom[1][0]], [rectgeom[0][1], rectgeom[1][1]]
        width = int(abs(rectgeom[1][0] - rectgeom[0][0]) / size)
        height = int(abs(rectgeom[1][1] - rectgeom[0][1]) / size)

        #Generating nodes
        for i in range(min(xlist), max(xlist) + 1, size):
            for j in range(min(ylist), max(ylist) + 1, size):
                nodes.append(self.prepclass.node(i, j))

        #Generating elements on nodes
        for i in range(width):
            for j in range(self.nodenum_4nextmesh, self.nodenum_4nextmesh + height):
                eles.append(self.prepclass.element(j + 1 + ((height + 1) * i),
                                                   j + 1 + ((height + 1) * (i + 1)),
                                                   j + ((height + 1) * (i + 1)),
                                                   j + ((height + 1) * i)))
        self.meshdict[self.meshnum] = [eles, 0, 0, size, nodes] #0 are placeholders for material
        self.nodenum_4nextmesh += 1 + nodes[-1] - nodes[0]
        return self.meshnum

    def mergenodes(self, nodenum_toleft, nodenum_tomerge):
        """Merges two nodes into one"""
        self.prepclass.swap_nodenum(nodenum_toleft, nodenum_tomerge)
        self.prepclass.deletenode(nodenum_tomerge)

    def assignproperty(self, meshnum, young, poisson, thickness):
        """Assigns material parameters and thickness to mesh"""
        self.meshdict[meshnum][1] = [young, poisson]
        self.meshdict[meshnum][2] = thickness

    def selectnode(self, x, y, meshnum=False):
        """Returns node number of given position and in given mesh"""
        nodenum_list = self.prepclass.nodes_bycoords(x, y)
        for nodenum in nodenum_list:
            if meshnum:
                if nodenum in self.meshdict[meshnum][4]:
                    return nodenum
            elif meshnum is False:
                return nodenum_list[0]
