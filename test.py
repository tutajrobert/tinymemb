import tinymemb

STEEL = [200e3, .3]

GEOM = tinymemb.Geom()
kp1 = (0, 0)
kp2 = (6, 8)
kp3 = (1, 0)
kp4 = (5, -6)

sup1 = GEOM.supele(kp1, kp2)
sup2 = GEOM.supele(kp3, kp4)

MESH = tinymemb.Mesh(GEOM)
mesh1 = MESH.generate(sup1, size=2)
mesh2 = MESH.generate(sup2, size=1)

MESH.nmerge(2, 0)
MESH.nmerge(4, 0)
MESH.assignprop(mesh1, STEEL, thickness=3)
MESH.assignprop(mesh2, STEEL, thickness=1)
#PREP.info()

SOLV = tinymemb.Solv(MESH)
SOLV.build()
SOLV.support(0, 8)
SOLV.support(2, 8)
SOLV.support(4, 8)
SOLV.support(6, 8)
SOLV.force(5, -6, 1000)

res = SOLV.solve()
print(min(res))
print(max(res))

"""
    a = [0, None, 3]
    b = [0, 5, 3]
    c = [None, None, None]
    d = [4, 4, 4]
    e = [0, None, 0]
     
    def cmp(a, b):
    	return all([x[0] == x[1] if x[0] is not None and x[1] is not None else True for x in zip(a,b)])
     
    print(cmp(a, b))
    print(cmp(b, b))
    print(cmp(b, e))
    print(cmp(b, c))
    print(cmp(d, e))
"""