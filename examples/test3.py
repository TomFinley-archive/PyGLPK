import glpk, random

def relaxed(similarities):
    def lin(a,b): return (a*(a+1))/2 + b
    size = len(similarities)
    lp = glpk.LPX()
    lp.params.msglev = 0
    lp.cols.add(size * (size+1) / 2)
    for col in lp.cols: col.bounds = 0,1
    for i in xrange(size):
        ind1 = lin(i,i)
        for j in xrange(i+1):
            ind = lin(i,j)
            lp.obj[ind] = similarities[i][j]
            if j==i: continue
            ind2 = lin(j,j)
            lp.rows.add(6)
            for row in lp.rows[-6:-3]: row.bounds = None,0
            lp.rows[-6].name = '%d,%d' %(i,j)
            lp.rows[-5].matrix = [(ind,1), (ind1,-1)]
            lp.rows[-4].matrix = [(ind,1), (ind2,-1)]
            for row in lp.rows[-3:]: row.bounds = None,1
            lp.rows[-3].matrix = [(ind,-1), (ind1,1), (ind2,1)]
            lp.rows[-2].matrix = [(ind,1), (ind1,-1), (ind2,1)]
            lp.rows[-1].matrix = [(ind,1), (ind1,1), (ind2,-1)]
    lp.obj.maximize = True
    lp.adv_basis()
    lp.simplex()

    label_matrix = [[col.value for col in lp.cols[lin(i,0):lin(i+1,0)]]
                    for i in xrange(size)]
    return label_matrix

def randomsim(size=5,randgen=random):
    """Create a random similarity matrix of given size."""
    sim = [[randgen.uniform(-1,1) for c2 in xrange(c+1)] for c in xrange(size)]
    return sim

random.seed(0)
for i in xrange(10000):
    s = randomsim(5)
    labs = relaxed(s)
