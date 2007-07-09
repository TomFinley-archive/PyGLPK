import glpk, random, operator, gc

random.seed(1)

# ---------------------------------------

print 'Testing adding rows and columns, and iterating on them...'

lp = glpk.LPX()
assert len(lp.rows)==lp.rows.add(3)==len(lp.rows)-3
assert len(lp.rows)==lp.rows.add(2)==len(lp.rows)-2
try:
    lp.cols.add(-2) # We want this to fail.
    assert False
except ValueError:
    pass
assert len(lp.rows)==5 and len(lp.cols)==0
assert len(lp.cols)==lp.cols.add(4)==len(lp.cols)-4
assert len(lp.cols)==4
assert len([r for r in lp.rows])==5 and len([c for c in lp.cols])==4

# ---------------------------------------

print 'Testing changing the objective function values...'

coefs = [7,2,9,-9]
lp.obj[:]=coefs
assert [lp.obj[c.index] for c in lp.cols]==coefs
lp.obj.shift = 3
assert lp.obj[None]==3
assert lp.obj[:2]==coefs[:2] and lp.obj[-3:]==coefs[-3:]
assert lp.obj[1,None,-1,2]==[2,3,-9,9]

# ---------------------------------------

print 'Testing naming objects...'

lp = glpk.LPX()

# Try setting and unsetting the linear program name.
assert lp.name == None
lp.name = 'albert'; assert lp.name == 'albert'
lp.name = 'bruce';  assert lp.name == 'bruce'
del lp.name;        assert lp.name == None

# Try setting and unsetting the objective function name.
assert lp.obj.name == None
lp.obj.name = 'nancy';  assert lp.obj.name == 'nancy'
lp.obj.name = 'olivia'; assert lp.obj.name == 'olivia'
del lp.obj.name;        assert lp.obj.name == None

lp.rows.add(60)
lp.cols.add(40)
# Ensure that all row and column names are None.
for r in lp.rows: assert r.name==None
for c in lp.cols: assert c.name==None
# Choose some random ones to name, and names.
rowi2name = dict([(r,chr(ord('A')+w)) for w,r in enumerate(
    random.sample(xrange(len(lp.rows)),26))])
coli2name = dict([(c,chr(ord('a')+w)) for w,c in enumerate(
    random.sample(xrange(len(lp.cols)),26))])
# Name them.
for r,n in rowi2name.items(): lp.rows[r].name = n
for c,n in coli2name.items(): lp.cols[c].name = n
# Go through all.  Check that those that are not named should not be,
# and those that are named have the right name.
for r in lp.rows: assert r.name==rowi2name.get(r.index,None)
for c in lp.cols: assert c.name==coli2name.get(c.index,None)
# Another related check.
assert dict([(r.index,r.name) for r in lp.rows if r.name])==rowi2name
assert dict([(c.index,c.name) for c in lp.cols if c.name])==coli2name
# Now test deleting all these names.
for r in lp.rows: del r.name; assert r.name==None
for c in lp.cols: del c.name; assert c.name==None

# ---------------------------------------

print 'Testing deletion of rows and columns...'

lp = glpk.LPX()
lp.rows.add(10)
lp.cols.add(10)
for r in lp.rows: r.name=chr(ord('a')+r.index)
for c in lp.cols: c.name=chr(ord('n')+c.index)
assert lp.rows[4].name=='e' and lp.rows[-2].name=='i' and lp.rows[8].name=='i'


# ---------------------------------------

print 'Testing setting up and querying the constraint matrix...'

# Make this work even in Python 2.3...
try: set
except NameError:
    import sets
    set = sets.Set

lp = glpk.LPX()
lp.rows.add(30)
lp.cols.add(50)

entries = [(r.index,c.index,random.uniform(-100,100)) # Create random entry
           for r in lp.rows for c in lp.cols          # for any row and col
           if random.random()>0.5]                    # with 50% probability.
entries = [e for e in entries if e[-1]]               # Ensure no zeros.

lp.matrix = entries
assert set(lp.matrix) == set(entries)
assert lp.nnz == len(entries)
rc2value = dict([((r,c),v) for r,c,v in entries])
lp.rows.add(3)

runningset = set()
for row in lp.rows:
    assert row.nnz == len([e for e in entries if e[0]==row.index])
    for cind, val in row.matrix:
        assert rc2value[row.index,cind] == val
        runningset.add((row.index,cind,val))
for col in lp.cols:
    assert col.nnz == len([e for e in entries if e[1]==col.index])
    for rind, val in col.matrix:
        assert rc2value[rind,col.index] == val
        runningset.add((rind,col.index,val))
assert set(entries)==runningset
lp.matrix = lp.matrix
assert set(lp.matrix) == set(entries)

# ---------------------------------------

print 'Testing the garbage collector...'

del lp
# Ensure we start with an empty garbage collector, and start saving references.
gc.collect()
gcd_flags = gc.get_debug()
gc.set_debug(gc.DEBUG_SAVEALL)
numcreates = 100
# Just create some meaningless problem.
for i in xrange(numcreates):
    lp = glpk.LPX()
    lp.name = 'some name'
    lp.rows.add(2)
    lp.cols.add(5)
    lp.rows[0].matrix = [1,2,3,4,5]
    lp.rows[1].name = 'second row'
    lp.obj[:] = [6,7,8,9,10]
    del lp
gc.collect()

type2count = {}
# Count each type of object that was garbage collected.
for ob in gc.garbage:
    type2count[type(ob)] = type2count.get(type(ob),0)+1
# Stop "debugging" now.
gc.set_debug(gcd_flags)
del gc.garbage[:]
gc.collect()

assert type2count.get(glpk.LPX,   0)==numcreates   # This many LPs created
assert type2count.get(glpk.BarCollection,0)==2*numcreates
                                                   # Twice as many bar collctns
assert type2count.get(glpk.Bar,   0)==0            # GC doesn't handle these
assert type2count.get(glpk.Objective,0)==numcreates# Each has an objective
assert type2count.get(glpk.Params,0)==numcreates   # and a params

# ---------------------------------------

print 'Testing simplex and interior point solver on random problem...'

# Approximate equality test.
def aeq(a, b): return abs(a-b)<1e-5
# Set the seed so we know which problem will be generated.
random.seed(1)
# Generate the problem.
lp = glpk.LPX()
lp.params.msglev = 0
lp.cols.add(10)
lp.rows.add(10)
lp.matrix = [random.uniform(-5,5) for i in xrange(len(lp.cols)*len(lp.rows))]
lp.obj[:] = [random.uniform(-5,5) for i in xrange(len(lp.cols))]
for r in lp.rows: r.bounds = -random.random(), random.random()
for c in lp.cols: c.bounds = -random.random(), random.random()
for solver in [lp.simplex, lp.interior]:
    lp.obj.maximize = False
    assert lp.status=='undef'
    assert solver()==None # Compare against known values.
    assert lp.status=='opt'
    assert aeq(lp.obj.value, -1.67146257641)
    assert aeq(sum([r.primal for r in lp.rows]), -0.27821176010)
    assert aeq(sum([r.dual   for r in lp.rows]),  1.73287645801)
    assert aeq(sum([r.primal for c in lp.cols]),  6.24802084152)
    assert aeq(sum([r.dual   for c in lp.cols]), -3.25513527209)
    lp.obj.maximize = True
    assert lp.status=='undef'
    assert solver()==None # Again compare against known values.
    assert lp.status=='opt'
    assert aeq(lp.obj.value, 1.57305544337)
    assert aeq(sum([r.primal for r in lp.rows]),  2.92436965935)
    assert aeq(sum([r.dual   for r in lp.rows]),  3.28775529757)
    assert aeq(sum([r.primal for c in lp.cols]), -1.09488627294)
    assert aeq(sum([r.dual   for c in lp.cols]), -2.96468593544)

# ---------------------------------------

print 'Testing simplex solver on simple 2D grid problem...'

lp = glpk.LPX()
lp.params.msglev = 0
lp.cols.add(2)
lp.obj[0,1] = 1
# Try very simple rules.
x1, x2 = lp.cols[0], lp.cols[1] # For convenience...
x1.name, x2.name = 'x1', 'x2'
x1.bounds = None, 1
x2.bounds = None, 2
lp.obj.maximize = True
assert lp.simplex()==None and lp.status=='opt'
assert aeq(lp.obj.value, 3) and aeq(x1.primal, 1) and aeq(x2.primal, 2)
# Now try pushing it into unbounded territory.
lp.obj[0] = -1
assert lp.simplex()==None and lp.status=='unbnd'
# Redefine the bounds so it is bounded again.
x1.bounds = -1, None
assert lp.simplex()==None and lp.status=='opt'
assert aeq(lp.obj.value, 3) and aeq(x1.primal, -1) and aeq(x2.primal, 2)
# Now add in a row constraint forcing it to a new point, (1/2, 2).
lp.rows.add(1)
lp.rows[0].matrix = [-2, 1]
lp.rows[0].bounds = None, 1
assert lp.simplex()==None and lp.status=='opt'
assert aeq(lp.obj.value, 1.5) and aeq(x1.primal, .5) and aeq(x2.primal, 2)
# Now add in another forcing it to point (1/4, 3/2).
lp.rows.add(1)
lp.rows[-1].matrix = [-2, -1]
lp.rows[-1].bounds = -2, None
assert lp.simplex()==None and lp.status=='opt'
assert aeq(lp.obj.value, 1.25) and aeq(x1.primal, .25) and  aeq(x2.primal, 1.5)
# Now go for the gusto.  Change a column constraint to force infeasibility.
x2.bounds = 2, None # Instead of x2<=2, must now be >=2!  Tee hee.
assert lp.simplex()==None and lp.status=='nofeas'
# By removing the first row constraint, we allow opt point (-1,4).
del lp.rows[0]
lp.std_basis()
assert lp.simplex()==None and lp.status=='opt'
assert aeq(lp.obj.value, 5) and aeq(x1.primal, -1) and aeq(x2.primal, 4)

# ---------------------------------------

print 'Testing MIP on solving SAT problems...'

def solve_sat(expression):
    if len(expression)==0: return [] # Trivial case.  Otherwise count vars.
    numvars = max([max([abs(v) for v in clause]) for clause in expression])
    lp = glpk.LPX()                  # Construct an empty linear program.
    lp.params.msglev = 0             # Stop the annoying output.
    lp.cols.add(2*numvars)           # As many columns as there are literals.
    for col in lp.cols:              # Literal must be between false and true.
        col.bounds = 0.0, 1.0
    def lit2col(lit):                # Function to compute column index.
        return [2*(-lit)-1,2*lit-2][lit>0]
    for i in xrange(1, numvars+1):   # Ensure "oppositeness" of literals.
        lp.rows.add(1)
        lp.rows[-1].matrix = [(lit2col(i), 1.0), (lit2col(-i), 1.0)]
        lp.rows[-1].bounds = 1.0     # Must sum to exactly 1.
    for clause in expression:        # Ensure "trueness" of each clause.
        lp.rows.add(1)
        lp.rows[-1].matrix = [(lit2col(lit), 1.0) for lit in clause]
        lp.rows[-1].bounds = 1, None # At least one literal must be true.
    retval = lp.simplex()            # Try to solve the relaxed problem.
    assert retval == None            # Should not fail in this fashion.
    if lp.status!='opt': return None # If no relaxed solution, no exact sol.

    lp.kind = int                    # Set kind to MIP problem
    for col in lp.cols: col.kind = int
    retval = lp.integer()            # Try to solve this integer problem.
    assert retval == None            # Should not fail in this fashion.
    if lp.status != 'opt': return None
    return [col.value > 0.99 for col in lp.cols[::2]]

def verify(expression, assignment):  # Test to check validity of assignment.
    return len(expression)==sum([reduce(max, [(abs(l)==l)==assignment[
        abs(l)-1] for l in c]) for c in expression]) # How ugly!!!

# Representations of CNFs.  Tuples are disjunctions.  <0 indicates negation.
# This one is possible.
e1=[(1,-3,-4),(-1,2,-4),(1,2,3),(1,-2,4),(1,2,-3),(-1,2,3),(1,3,4),
    (-1,-2,3),(1,3,-4),(-2,-3,-4),(-2,3,-4),(2,3,4)]
# This one is impossible.
e2=[(1,-2,3),(-1,2,4),(-1,-2,-4),(1,2,4),(1,-3,4),(1,-2,-4),(1,2,3),
    (1,3,-4),(-2,-3,-4),(2,-3,-4),(-2,3,-4),(-2,-3,4),(-2,3,4),(2,3,-4)]
assert verify(e1, solve_sat(e1))     # Should have a valid solution.
assert solve_sat(e2)==None           # Should have no solution.

# ---------------------------------------
if glpk.version < (4, 7):
    print 'NOTE: name indexing not supported for GLPK %d.%d,'%glpk.version,
    print 'skipping those tests'
else:
    print 'Testing string name indices for row and column access...'

    lp = glpk.LPX()

    lp.cols.add(5)

    lp.cols[1].name = 'x1'
    lp.cols[4].name = 'x2'

    assert (lp.cols['x1'].index, lp.cols['x2'].index) == (1,4)
    assert [c.index for c in lp.cols['x1','x2','x2']] == [1,4,4]
    del lp.cols[2:4]
    assert [c.index for c in lp.cols['x1','x2','x2']] == [1,2,2]

    del lp.cols['x2','x1']
    try:
        print lp.cols['x2']
        assert False
    except KeyError:
        pass

    lp.rows.add(4)

    lp.rows[2].name = 'bob'
    assert lp.rows['bob'].index == 2
    lp.rows['bob'].name = 'stan'
    try:
        print lp.rows['bob']
        assert False
    except KeyError:
        pass
    assert 'bob' not in lp.rows
    assert 'stan' in lp.rows

# ---------------------------------------

    print 'Testing string name indices for matrix and obj function setting...'

    lp = glpk.LPX()

    lp.cols.add(4)
    for c in lp.cols: c.name = chr(ord('a')+c.index)
    lp.rows.add(3)
    for r in lp.rows: r.name = chr(ord('i')+r.index)

    lp.rows['j'].matrix = [('b',.2), .4, ('a',.6)]
    assert set(lp.matrix)==set([(1,1,.2), (1,2,.4), (1,0,.6)])
    del lp.matrix
    assert lp.nnz==0 and len(lp.matrix)==0
    lp.matrix = [(2,'a',.2), ('i',-1,.4), ('j','c',.6)]
    assert set(lp.matrix)==set([(2,0,.2), (0,3,.4), (1,2,.6)])
    lp.obj[:] = [-4,-5,-6,-7]

    assert lp.obj['b']==-5
    assert lp.obj['b','a',1,'c',-1]==[-5,-4,-5,-6,-7]
    assert lp.obj[::-1]==[-7,-6,-5,-4]

    lp.obj['b']=5
    assert lp.obj[1]==5
    lp.obj['a','c',None,-1] = 1,2,3,4
    assert lp.obj[0]==1 and lp.obj[2]==2 and lp.obj[None]==3 and lp.obj[3]==4
    assert lp.obj['a','c',None,-1]==[1,2,3,4]
    assert lp.obj.shift==3

# ---------------------------------------

print 'Tests succeeded!'
