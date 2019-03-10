import tinymemb

STEEL_modulus = 200e3
STEEL_poissratio = .3

GEOM = tinymemb.Geom()
p1 = GEOM.point(0, 0)
p2 = GEOM.point(6, 8)
p3 = GEOM.point(1, 0)
p4 = GEOM.point(5, -6)
rec1 = GEOM.rectangle(p1, p2)
rec2 = GEOM.rectangle(p3, p4)

MESH = tinymemb.Mesh(GEOM)
mesh1 = MESH.generate(rec1, size=2)
mesh2 = MESH.generate(rec2, size=1)

node1 = MESH.selectnode(2, 0, mesh1)
node2 = MESH.selectnode(2, 0, mesh2)
node3 = MESH.selectnode(4, 0, mesh1)
node4 = MESH.selectnode(4, 0, mesh2)

#MESH.mergenodes(node1, node2)
#MESH.mergenodes(node3, node4)
MESH.assignproperty(mesh1, young=200e3, poisson=.3, thickness=3)
MESH.assignproperty(mesh2, young=200e3, poisson=.3, thickness=1)

SOLV = tinymemb.Solv(MESH)
SOLV.build()
SOLV.connector(node1, node2)
#SOLV.connector(node1, node4)
SOLV.connector(node3, node4)
#SOLV.connector(node2, node3)
SOLV.support(0, 8)
SOLV.support(2, 8)
SOLV.support(4, 8)
SOLV.support(6, 8)
SOLV.force(5, -6, 1000)

res = SOLV.solve()

SOLV2 = tinymemb.Solv(MESH)
SOLV2.build()
#SOLV2.connector(node1, node2)
#SOLV.connector(node1, node4)
#SOLV2.connector(node3, node4)
force1 = SOLV2.backsolve4force(2, 0, res)
force2 = SOLV2.backsolve4force(4, 0, res)
force3 = SOLV2.backsolve4force(5, -6, res)
react = SOLV2.backsolve4force(0, 8, res)
react2 = SOLV2.backsolve4force(2, 8, res)
react3 = SOLV2.backsolve4force(4, 8, res)
react4 = SOLV2.backsolve4force(6, 8, res)

print(min(res))
print(max(res))
print(force1)
print(force2)
print(force3)
print(react, react2, react3, react4)
print(react + react2 + react3 + react4)

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