import tinymemb

GEOM = tinymemb.Geom()
p1 = GEOM.point(0, 120)
p4 = GEOM.point(40, 0)
p6 = GEOM.point(160, 40)

rec1 = GEOM.rectangle(p1, p4)
rec2 = GEOM.rectangle(p4, p6)

sz = 20

MESH = tinymemb.Mesh(GEOM)
mesh1 = MESH.generate(rec1, size=sz)
mesh2 = MESH.generate(rec2, size=sz)

#line select
nodes_mesh1 = []
nodes_mesh2 = []
for i in range(0, 41, sz):
    nodes_mesh1.append(MESH.selectnode(40, i, mesh1))
    nodes_mesh2.append(MESH.selectnode(40, i, mesh2))
"""
for node1, node2 in zip(nodes_mesh1, nodes_mesh2):
    MESH.mergenodes(node1, node2)
"""
MESH.assignproperty(mesh1, young=200e3, poisson=.3, thickness=2)
MESH.assignproperty(mesh2, young=200e3, poisson=.3, thickness=2)

SOLV = tinymemb.Solv(MESH)
SOLV.build()

for node1, node2 in zip(nodes_mesh1, nodes_mesh2):
    SOLV.connector(node1, node2, k=1e9)

#SOLV.print_matrix()

nodes_support = []
for i in range(0, 41, sz):
    nodes_support.append(MESH.selectnode(160, i))
    SOLV.support(nodes_support[-1])

node_f = MESH.selectnode(40, 120)
SOLV.force(node_f, xforce=0, yforce=1e3)

results = SOLV.solve()

"""
reaction_x, reaction_y = [], []
for node in nodes_support:
    print(SOLV.backsolve_4force(node))
    reaction_x.append(SOLV.backsolve_4force(node)[0])
    reaction_y.append(SOLV.backsolve_4force(node)[1])

print(SOLV.backsolve_4force(node_f))

print(sum(reaction_x))
print(sum(reaction_y))
"""

POST = tinymemb.Post(SOLV)
res1 = POST.disp(MESH.selectnode(0, 40))
res2 = POST.disp(MESH.selectnode(40, 40))


print(res1)
print(res2)

print(SOLV.conn_force(2))
print(SOLV.conn_force(1))
print(SOLV.conn_force(0))