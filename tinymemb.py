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
    def __init__(self, meshclass):
        self.prepclass = meshclass.prepclass
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
        nnum = self.prepclass.nodes_bycoords(x, y)[0]
        dofs = [self.nod2dofmap[nnum], self.nod2dofmap[nnum] + 1]
        for dof in dofs:
            self.gklist[dof] = [0 for i in range(self.nodesnum * 2)]
            for row in self.gklist:
                row[dof] = 0
            self.gklist[dof][dof] = 1

    def force(self, x, y, force_value):
        nnum = self.prepclass.nodes_bycoords(x, y)[0]
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
        nnum = self.prepclass.nodes_bycoords(x, y)[0]
        print(nnum)
        dof = self.nod2dofmap[nnum]
        force = 0
        for i in range(0, self.nodesnum * 2):
            if self.gklist[dof][i] * dlist[i]:
                force += self.gklist[dof][i] * dlist[i]
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
        
    def selectnode(self, x, y, meshnum):
        """Returns node number of given position and in given mesh"""
        nodenum_list = self.prepclass.nodes_bycoords(x, y)
        for nodenum in nodenum_list:
            if nodenum in self.meshdict[meshnum][4]:
                return nodenum
