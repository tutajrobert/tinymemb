import tinymemb

GEOM = tinymemb.Geom()
p1 = GEOM.point(1, 0)
p2 = GEOM.point(5, 6)
p3 = GEOM.point(1, 0)
p4 = GEOM.point(5, -6)
rec1 = GEOM.rectangle(p1, p2)
rec2 = GEOM.rectangle(p3, p4)

MESH = tinymemb.Mesh(GEOM)
mesh1 = MESH.generate(rec1, size=1)
mesh2 = MESH.generate(rec2, size=1)

node1 = MESH.selectnode(2, 0, mesh1)
node2 = MESH.selectnode(2, 0, mesh2)
node3 = MESH.selectnode(4, 0, mesh1)
node4 = MESH.selectnode(4, 0, mesh2)

#MESH.mergenodes(node1, node2)
#MESH.mergenodes(node3, node4)
MESH.assignproperty(mesh1, young=200e3, poisson=.5, thickness=2)
MESH.assignproperty(mesh2, young=200e3, poisson=.5, thickness=2)

SOLV = tinymemb.Solv(MESH)
SOLV.build()
SOLV.connector(node1, node2)
SOLV.connector(node3, node4)
nodesup2 = MESH.selectnode(1, 6)
nodesup3 = MESH.selectnode(2, 6)
nodesup4 = MESH.selectnode(3, 6)
nodesup5 = MESH.selectnode(4, 6)
nodesup6 = MESH.selectnode(5, 6)
nodesup_list = [nodesup2, nodesup3, nodesup4, nodesup5, nodesup6]
for nodesup in nodesup_list:
    SOLV.support(nodesup)
nodeforce1 = MESH.selectnode(5, -6)
nodeforce2 = MESH.selectnode(1, -6)
nodeforce3 = MESH.selectnode(3, -6)
#SOLV.force(nodeforce1, xforce=1000, yforce=0)
#SOLV.force(nodeforce2, xforce=1000, yforce=0)
SOLV.force(nodeforce3, xforce=1000, yforce=0)

results = SOLV.solve()

force1 = SOLV.backsolve_4force(node2)
force2 = SOLV.backsolve_4force(node4)
force1b = SOLV.backsolve_4force(node1)
force2b = SOLV.backsolve_4force(node3)
force3 = SOLV.backsolve_4force(nodeforce1)
react2 = SOLV.backsolve_4force(nodesup2)
react3 = SOLV.backsolve_4force(nodesup3)
react4 = SOLV.backsolve_4force(nodesup4)

POST = tinymemb.Post(SOLV)

res1 = POST.disp(nodeforce1)
res2 = POST.disp(nodeforce2)

#print(res1)
#print(res2)

print(force1)
print(force1b)
print(force2)
print(force2b)
#print(force3)
#print(react, react2, react3, react4)