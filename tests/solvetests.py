"""Tests for the solver itself."""

from testutils import *

class TwoDimensionalTest(unittest.TestCase):
    def setUp(self):
        self.lp = LPX()
        self.lp.params.msglev = 0

    def runTest(self):
        """Test repeatedly solving an LP with evolving constraints."""
        lp = self.lp
        # Set up the rules of the problem.
        lp.cols.add(2)
        lp.obj[0,1] = 1
        # Try very simple rules.
        x1, x2 = lp.cols[0], lp.cols[1] # For convenience...
        x1.name, x2.name = 'x1', 'x2'
        x1.bounds = None, 1
        x2.bounds = None, 2
        lp.obj.maximize = True
        self.assertEqual(None, lp.simplex())
        self.assertEqual('opt', lp.status)
        self.assertAlmostEqual(lp.obj.value, 3)
        self.assertAlmostEqual(x1.primal, 1)
        self.assertAlmostEqual(x2.primal, 2)
        # Now try pushing it into unbounded territory.
        lp.obj[0] = -1
        self.assertEqual(None, lp.simplex()) # No error?
        self.assertEqual('unbnd', lp.status)
        # Redefine the bounds so it is bounded again.
        x1.bounds = -1, None
        self.assertEqual(None, lp.simplex())
        self.assertEqual('opt', lp.status)

        self.assertAlmostEqual(lp.obj.value, 3)
        self.assertAlmostEqual(x1.primal, -1)
        self.assertAlmostEqual(x2.primal, 2)
        # Now add in a row constraint forcing it to a new point, (1/2, 2).
        lp.rows.add(1)
        lp.rows[0].matrix = [-2, 1]
        lp.rows[0].bounds = None, 1
        self.assertEqual(None, lp.simplex())
        self.assertEqual('opt', lp.status)
        self.assertAlmostEqual(lp.obj.value, 1.5)
        self.assertAlmostEqual(x1.primal, .5)
        self.assertAlmostEqual(x2.primal, 2)
        # Now add in another forcing it to point (1/4, 3/2).
        lp.rows.add(1)
        lp.rows[-1].matrix = [-2, -1]
        lp.rows[-1].bounds = -2, None

        self.assertEqual(None, lp.simplex())
        self.assertEqual('opt', lp.status)

        self.assertAlmostEqual(lp.obj.value, 1.25)
        self.assertAlmostEqual(x1.primal, .25)
        self.assertAlmostEqual(x2.primal, 1.5)
        # Now go for the gusto.  Change column constraint, force infeasibility.
        x2.bounds = 2, None # Instead of x2<=2, must now be >=2!  Tee hee.
        self.assertEqual(None, lp.simplex())
        self.assertEqual('nofeas', lp.status)
        # By removing the first row constraint, we allow opt point (-1,4).
        del lp.rows[0]
        lp.std_basis()

        self.assertEqual(None, lp.simplex())
        self.assertEqual('opt', lp.status)

        self.assertAlmostEqual(lp.obj.value, 5)
        self.assertAlmostEqual(x1.primal, -1)
        self.assertAlmostEqual(x2.primal, 4)

class SatisfiabilityMIPTest(unittest.TestCase):
    @classmethod
    def solve_sat(self, expression):
        """Attempts to satisfy a formula of conjunction of disjunctions.

        If there are n variables in the expression, this will return a
        list of length n, all elements booleans.  The truth of element i-1
        corresponds to the truth of variable i.

        If no satisfying assignment could be found, None is returned."""
        # Trivial boundary case.
        if len(expression)==0: return []
        # Now many variables do we have?
        numvars = max(max(abs(v) for v in clause) for clause in expression)

        # Construct an empty linear program.
        lp = LPX()

        # The output GLPK produces is rather annoying.
        lp.params.msglev = 0

        # We want twice as many columns (LP variables) as there are
        # logical variables in the expression: one column for each
        # positive literal, and one column for each negative literal.  A
        # literal is "true" if its column holds 1, false if 0.
        lp.cols.add(2*numvars)

        # Bound all columns (LP variables) to have value between 0 and 1.
        for col in lp.cols:
            col.bounds = 0.0, 1.0

        # Let us suppose that literals x_i and !x_i correspond to columns
        # 2*i-2 and 2*i-1 respectively.  (So, columns 0 and 1 correspond
        # to the truth of x_1 and !x_1, respectively, columns 2 and 3 to
        # the truth of x_2 and !x_2, etc.)

        # Let's just define a helper function that will perform the
        # mapping of literal identifier to column number.  (So 1 maps to
        # 0, -1 to 1, 2 to 2, -2 to 3, 3 to 4, etc.)
        def lit2col(lit):
            return 2*lit-2 if lit>0 else 2*(-lit)-1

        # Here we do two things: we set the names of the column variables
        # appropriately (not required, but fun), and define a row in a
        # constraint matrix indicating that a literal and its complement
        # must sum to exactly 1 since exactly one of the two literals must
        # be true.
        for i in xrange(1, numvars+1):
            lp.cols[lit2col( i)].name =  'x_%d'%i
            lp.cols[lit2col(-i)].name = '!x_%d'%i
            lp.rows.add(1)
            lp.rows[-1].matrix = [(lit2col(i), 1.0), (lit2col(-i), 1.0)]
            lp.rows[-1].bounds = 1.0

        # We want to find an assignment of such variables so that each
        # disjunction (or-ing) is true.  A disjunction is true if one of
        # its literals is true.  Literal truth corresponds to 1 in the
        # corresponding column, so we can represent truth of a disjunction
        # by defining a constraint saying the sum of its literals must be
        # at least 1.  We do this for all clauses.
        for clause in expression:
            lp.rows.add(1)
            lp.rows[-1].matrix = [(lit2col(lit), 1.0) for lit in clause]
            lp.rows[-1].bounds = 1, None

        # Now we have the LP built.  Run the simplex algorithm.
        retval = lp.simplex()

        # If our iteration terminated prematurely, or if we do not have an
        # optimal solution, assume we have failed.
        if retval != None: return None
        if lp.status != 'opt': return None

        # Now switch this from a linear continuous problem to a
        # mixed-integer problem, and say we want all column variables
        # integer as well (that is, not just between 0 and 1, but exactly
        # 0 and 1).
        lp.kind = int
        for col in lp.cols:
            col.kind = int

        # Attempt to solve this MIP problem with the MIP solver.
        retval = lp.integer()

        # Again, only returns non-None on failure.
        if retval != None: return None
        if lp.status != 'opt': return None

        # We want to return a list of boolean quantities, where the first
        # element is True iff x_1 is true, the second element is True iff
        # x_2 is true, and so on.  Variable truth corresponds to the
        # literals represented in the even columns.  So, we can just go
        # over all of the even columns (every 2 counting from 0), see if
        # the value is close to 1, and use that as our variable
        # assignment.
        return [col.value > 0.99 for col in lp.cols[::2]]

    @classmethod
    def verify(self, expression, assignment):
        """Get the truth of an expression given a variable truth assignment.

        This will return true only if this is a satisfying assignment."""
        # Each clause must be true.
        for clause in expression:
            # For a disjunctive clause to be true, at least one of its
            # literals must be true.
            lits_true = 0
            for lit in clause:
                if (lit>0) == (assignment[abs(lit)-1]):
                    lits_true += 1
            if lits_true == 0:
                return False
        return True

    def testSolvableSat(self):
        """Test that a solvable satisfiability problem can be solved."""
        exp = [(3, -2, 5), (1, 2, 5), (1, 3, -4), (-1, -5, -3), (-4, 2, 1),
               (5, 2, -1), (5, 2, 1), (-4, -1, -5), (4, 5, -1), (3, 1, -2),
               (1, 5, -3), (-5, -3, -1), (-4, -5, -3), (-3, -5, -2),
               (4, -5, 3), (1, -2, -5), (1, 4, -3), (4, -1, -3),
               (-5, -1, 3), (-2, -4, -5)]
        solution = self.solve_sat(exp)
        # First assert that this is a solution.
        self.assertNotEqual(solution, None)
        self.failUnless(self.verify(exp, solution))

    def testInsolvableSat(self):
        """Test that an unsolvable satisfiability problem can't be solved."""
        exp = [(1,2), (1,-2), (-1,2), (-1,-2)]
        solution = self.solve_sat(exp)
        self.assertEqual(solution, None)
