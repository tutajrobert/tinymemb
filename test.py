import tinymemb

PREP = tinymemb.Prep()

GEOM = tinymemb.Geom()
kp1 = GEOM.kpoint(0, 0)
kp2 = GEOM.kpoint(2, 4)
kp3 = GEOM.kpoint(2, -4)
sup1 = GEOM.supele(kp1, kp2)
sup2 = GEOM.supele(kp1, kp3)

elast_mod = 205e3
mesh1 = PREP.mesh(sup1, 1, elast_mod)
mesh2 = PREP.mesh(sup2, 1, elast_mod)

PREP.nmerge(0, 0)
PREP.nmerge(2, 0)
#PREP.info()

SOLV = tinymemb.Solv(PREP)

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