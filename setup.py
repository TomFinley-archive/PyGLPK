from distutils.core import setup, Extension
import sys, os, os.path

source_roots = 'glpk lp barcol bar obj util params kkt'.split()

libdirs, incdirs = [], []

# Attempt to automatically detect which version of GLPK we are using.
try:
    fi, fo = os.popen4('glpsol -v')
    fi.close()
    tokens = fo.read().split()
    version_string = tokens[tokens.index('Version')+1]
    print 'Version GLPK %s detected!!' % version_string
    version = version_string.split('.')
    version = tuple([int(i) for i in version[:2]])

    # Try to get which is the executable path, and infer additional
    # library and include directories from there.
    glpsol_path = os.popen('which glpsol').read().strip()
    glpsol_path = os.path.abspath(glpsol_path)

    head, tail = os.path.split(glpsol_path)
    head, tail = os.path.split(head)
    libdirs.append(os.path.join(head, 'lib'))
    incdirs.append(os.path.join(head, 'include'))
    
except Exception, e:
    print 'Failed to run glpsol and extract version number!'.upper()
    print 'Just assuming GLPK 4.4.  (Note, some functions will not work.)'
    version = (4, 4)

# CUSTOM USER SECTION

# If necessary, manually change this tuple to reflect the actual
# version of GLPK, e.g., version x.y should have the tuple (x, y).
#version = (4, 15)

# If necessary, uncomment and set your library search path manually.
#libdirs = ['/some/path/to/lib']

# If necessary, uncomment and set your header search path manually.
#incdirs = ['/some/path/to/include']

# If necessary, uncomment and define what libraries to include.
#libs = []

# If necessary, uncomment and define extra object files to include.
extraobs = []

# END OF CUSTOM USER SECTION

if version < (4, 4):
    # Many key functions were renamed in 4.4 and later.
    print 'GLPK version appears to be earlier than 4.4!'
    print 'This probably won\'t work.'

# If the user did not define libraries themselves, set that up.  Only
# need to have GMP if we have a version that uses exact arithmetic.
try:
    libs
except NameError:
    # The user did not set that up.
    libs = ['glpk']
    if version >= (4, 13):
        libs.append('gmp')

incdirs.append('src')

# Define the version macros so that the C compiler can make decisions.
versions = [('GLPK_MAJOR_VERSION', str(version[0])),
            ('GLPK_MINOR_VERSION', str(version[1]))]
# Now, finally, define that module!
module1 = Extension(
    'glpk',
    sources = [os.path.join('src',r+'.c') for r in source_roots],
    define_macros = versions,
    library_dirs = libdirs, include_dirs = incdirs,
    libraries = libs, extra_objects = extraobs)

ld = """The PyGLPK module gives one access to the functionality
of the GNU Linear Programming Kit.  
"""

setup(name = 'glpk',
      version = '0.2',
      description = 'PyGLPK, a Python module encapsulating GLPK.',
      long_description = ld,
      author = 'Thomas Finley',
      author_email = 'tomf@cs.cornell.edu',
      url = 'http://www.cs.cornell.edu/~tomf/pyglpk/',
      license = 'GPL',
      classifiers = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Programming Language :: C',
    'Programming Language :: Python',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Software Development :: Libraries :: Python Modules' ],
      ext_modules = [module1])
