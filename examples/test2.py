import glpk            # Import the GLPK module

lp = glpk.LPX()        # Create empty problem instance
lp.name = 'sample'     # Assign symbolic name to problem
lp.obj.maximize = True # Set this as a maximization problem
lp.rows.add(3)         # Append three rows to this instance
for r in lp.rows:      # Iterate over all rows
	r.name = chr(ord('p')+r.index) # Name them p, q, and r
lp.rows[0].bounds = None, 100.0  # Set bound -inf < p <= 100
lp.rows[1].bounds = None, 600.0  # Set bound -inf < q <= 600
lp.rows[2].bounds = None, 300.0  # Set bound -inf < r <= 300
lp.cols.add(3)         # Append three columns to this instance
for c in lp.cols:      # Iterate over all columns
	c.name = 'x%d' % c.index # Name them x0, x1, and x2
	c.bounds = 0.0, None     # Set bound 0 <= xi < inf
lp.obj[:] = [ 10.0, 6.0, 4.0 ]   # Set objective coefficients
lp.matrix = [ 1.0, 1.0, 1.0,     # Set nonzero entries of the
             10.0, 4.0, 5.0,     #   constraint matrix.  (In this
              2.0, 2.0, 6.0 ]    #   case, all are non-zero.)
lp.simplex()           # Solve this LP with the simplex method
print 'Z = %g;' % lp.obj.value,  # Retrieve and print obj func value
print '; '.join(['%s = %g' % (c.name, c.primal) for c in lp.cols])
                       # Print struct variable names and primal values
