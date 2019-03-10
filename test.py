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
SOLV.connector(node3, node4)
nodesup1 = MESH.selectnode(0, 8)
nodesup2 = MESH.selectnode(2, 8)
nodesup3 = MESH.selectnode(4, 8)
nodesup4 = MESH.selectnode(6, 8)
nodesup_list = [nodesup1, nodesup2, nodesup3, nodesup4]
for nodesup in nodesup_list:
    SOLV.support(nodesup)
nodeforce = MESH.selectnode(5, -6)
SOLV.force(nodeforce, xforce=-1000, yforce=0)

results = SOLV.solve()

force1 = SOLV.backsolve_4force(node2)
force2 = SOLV.backsolve_4force(node4)
force3 = SOLV.backsolve_4force(nodeforce)
react = SOLV.backsolve_4force(nodesup1)
react2 = SOLV.backsolve_4force(nodesup2)
react3 = SOLV.backsolve_4force(nodesup3)
react4 = SOLV.backsolve_4force(nodesup4)

print(min(results))
print(max(results))
print(force1)
print(force2)
print(force3)
print(react, react2, react3, react4)
print(react + react2 + react3 + react4)