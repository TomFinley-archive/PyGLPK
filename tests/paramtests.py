"""Tests for setting parameters."""

from testutils import *
import random

class ParamsTestCase(Runner, unittest.TestCase):
    """Tests for setting and getting parameters of the LP."""
    def setUp(self):
        self.lp = LPX()
        param_names = [n for n in dir(self.lp.params) if n[0]!='_']
        param_names.remove('reset')
        self.param_name2default = dict([
            (n, getattr(self.lp.params, n)) for n in param_names])

        # INTEGER PARAMETER VALUES

        # Have nice lists of, for each parameter, valid values, and
        # values illegal because they are of bad type or value.
        p2v, p2bt, p2bv = {}, {}, {}
        # These have legal values 0,1
        for n in 'price'.split():
            p2v[n] = 0,1
            p2bv[n] = -100, -1, 2, 3, 4, 5
        # These have legal values 0,1,2
        for n in 'branch btrack mpsobj'.split():
            p2v[n] = 0,1,2
        # These have legal values 0,1,2,3
        for n in 'msglev scale'.split():
            p2v[n] = 0,1,2,3
        # These all have the following illegal typed values.
        illegals = 'hi!', [1], {2:5}, 0.5, complex(1,0), self.lp
        for n in p2v:
            p2bt[n] = illegals
        # Set their illegal values.
        illegal_vals = -100, -1, 2, 3, 4, 5, 1000
        for n in p2v:
            legals = set(p2v[n])
            p2bv[n] = tuple([v for v in illegal_vals if v not in legals])
        # These may be any positive integer.
        for n in 'outfrq'.split():
            p2v[n] = 1, 2, 3, 100, 1600, 80000, 4000000, 2000000000
            p2bt[n] = illegals
            p2bv[n] = 0, -1, -2, -3, -100, -100000
        # These may be any integer.
        for n in 'itlim'.split():
            p2v[n] = illegal_vals
            p2bt[n] = illegals
        # Integer in range 255.
        for n in 'usecuts'.split():
            p2v[n] = xrange(0, 0x100)
            p2bt[n] = illegals
            p2bv[n] = tuple(xrange(-10,0)) + tuple(xrange(0x100, 0x110))

        # FLOAT PARAMETER VALUES

        # Previous illegal types included floats.
        illegals = tuple([iv for iv in illegals if type(iv)!=float])
        # Floats without bound.
        for n in 'objll objul outdly tmlim'.split():
            p2v[n] = -1e100, -5.23e20, -2.231, -1, 0, 0.5, 1, 3.14159, 1e100
            p2bt[n] = illegals
        # Bounded in [0,1] range.
        for n in 'relax'.split():
            p2v[n] = 0, .2, .45678, .901, 1
            p2bv[n] = -5e6, -1, -0.01, 1.01, 1.5, 1000
            p2bt[n] = illegals
        # Bounded in double-epsilon to 0.001 range.
        for n in 'tolbnd toldj tolint tolobj tolpiv'.split():
            p2v[n] = 2.3e-16, 1e-8, 0.000523, 0.001
            p2bv[n] = -1000, -1e-8, 0, 0.0011, 0.1, 5.123, 1e100
            p2bt[n] = illegals

        # BOOLEAN PARAMETER VALUES

        # These have boolean values.  They have no illegal values
        # owing to PyGLPK's pecularity of having boolean values set to
        # the bool(value) of value, which is always defined.
        bname = 'dual mpsfree mpsinfo mpsorig mpsskip mpswide presol round'
        for n in bname.split():
            p2v[n] = False, True

        # Set up the mapping from parameter names to tuples of legal
        # values, and illegal params because of value and type.
        self.param_name2values = p2v
        self.param_name2bad_values = p2bv
        self.param_name2bad_type_values = p2bt

    def testIterationCountReadOnly(self):
        """The itcnt parameter should be read only."""
        # TypeError was raised in pre-2.5 versions of Python.
        # AttributeError is raised in 2.5-post versions of Python.
        self.assertRaises(Exception, self.runner,
                          'self.lp.params.itcnt=5')
        self.assertRaises(Exception, self.runner,
                          'self.lp.params.itcnt=-300')

    def testDefaultSetIsNotIllegal(self):
        """Sets all parameters to their default values."""
        p2dv = dict(self.param_name2default)
        del p2dv['itcnt']
        for pname, pvalue in p2dv.iteritems():
            setattr(self.lp.params, pname, pvalue)

    def listOfParameterPairs(self, pname2pvalues):
        """Reduce a parameter name, value-sequence mapping to a list.

        Given a mapping from parameter names to a parameter value
        sequence, construct a list of all parameter names and
        parameter values."""
        pp = [(name, value) for name, values in pname2pvalues.iteritems()
              for value in values]
        return pp

    def testSetValues(self):
        """Set many parameter values."""
        pp = self.listOfParameterPairs(self.param_name2values)
        random.Random(2).shuffle(pp)
        for n, v in pp:
            setattr(self.lp.params, n, v)
            self.assertEqual(getattr(self.lp.params, n), v)

    def testSetValuesThenReset(self):
        """Set many parameter values, then use reset() to set defaults."""
        pp = self.listOfParameterPairs(self.param_name2values)
        pp = [(n,v) for n,v in pp if self.param_name2default[n]!=v]
        random.Random(3).shuffle(pp)
        for n, v in pp:
            setattr(self.lp.params, n, v)
        self.lp.params.reset()
        for n, v in sorted(self.param_name2default.iteritems()):
            self.assertEqual(getattr(self.lp.params, n), v)

    def testSetBadValues(self):
        """Set many parameters values of bad value."""
        pp = self.listOfParameterPairs(self.param_name2bad_values)
        random.Random(4).shuffle(pp)
        for n, v in pp:
            self.assertRaises(ValueError, setattr, self.lp.params, n, v)

    def testSetBadTypeValues(self):
        """Set many parameters values of bad type."""
        pp = self.listOfParameterPairs(self.param_name2bad_type_values)
        random.Random(5).shuffle(pp)
        for n, v in pp:
            self.assertRaises(TypeError, setattr, self.lp.params, n, v)

class IntegerProgramTestCase(Runner, unittest.TestCase):
    """Tests for setting the program type (continuous/int)."""
    def setUp(self):
        self.lp = LPX()

    def testDefaultKind(self):
        """Test that the default kind is float (continuous)."""
        self.assertEqual(self.lp.kind, float)
        self.lp.cols.add(5)
        self.assertEqual([c.kind for c in self.lp.cols], [float]*5)
        self.lp.rows.add(4)
        self.assertEqual([r.kind for r in self.lp.rows], [float]*4)

    def testSettingKind(self):
        """Test setting the kind of the problem."""
        self.lp.kind = int
        self.assertEqual(self.lp.kind, int)
        self.lp.kind = float
        self.assertEqual(self.lp.kind, float)

    def testSettingColumnsIntBeforeProblemIsInt(self):
        """Test setting columns as ints before setting the problem int."""
        self.lp.cols.add(5)
        self.assertRaises(ValueError, self.runner,
                          'self.lp.cols[2].kind = int')

    def testSettingColumnsInt(self):
        """Test setting columns as int."""
        self.lp.cols.add(5)
        self.lp.kind = int
        for c in [-1, 1, 2]:
            self.lp.cols[c].kind = int
        self.assertEqual([c.kind for c in self.lp.cols],
                         [float, int, int, float, int])

    def testIntAndBinaryCounters(self):
        """Test nint and nbin counters."""
        self.lp.cols.add(5)
        self.lp.kind = int
        # Convenience function to check the values of nint and nbin.
        def nn(ni,nb): self.assertEqual((self.lp.nint, self.lp.nbin), (ni,nb))
        # Do a series of ops and watch how nint and nbin changes.
        nn(0,0)
        self.lp.cols[ 0].kind = int;    nn(1,0)
        self.lp.cols[ 4].kind = int;    nn(2,0)
        self.lp.cols[ 4].bounds = 0,2;  nn(2,0)
        self.lp.cols[-1].kind = int;    nn(2,0)
        self.lp.cols[ 0].kind = float;  nn(1,0)
        self.lp.cols[ 3].kind = int;    nn(2,0)
        self.lp.cols[ 4].bounds = 0,1;  nn(2,1)
        self.lp.cols[ 1].kind = int;    nn(3,1)
        self.lp.cols[ 3].bounds = 0,1;  nn(3,2)
        self.lp.cols[ 4].bounds = None; nn(3,1)
        del self.lp.cols[3];            nn(2,0)
        self.lp.cols[ 2].bounds = 0,1;  nn(2,0)
        self.lp.cols[ 2].kind = int;    nn(3,1)

    def testSettingsRowsInt(self):
        """Test setting rows as int, which should always fail."""
        self.lp.rows.add(5)
        self.assertRaises(ValueError, self.runner,
                          'self.lp.rows[2].kind = int')
        self.lp.kind = int
        self.assertRaises(ValueError, self.runner,
                          'self.lp.rows[4].kind = int')

    def testSettingsRowsWithBadObject(self):
        """Test setting columns as something neither int or float."""
        self.lp.cols.add(5)
        for kind in ['complex', 'bool', '"belay"', '2', '{1:3}']:
            self.assertRaises(ValueError, self.runner,
                              'self.lp.cols[2].kind = %r' % kind)

