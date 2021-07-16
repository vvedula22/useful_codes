# Header
import sympy as sp
sp.init_printing(use_unicode = True)

class bcolors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    ULINE   = '\033[4m'

# Define symbols
x = sp.symbols('x')

s3x = sp.sin(x) * sp.sin(x) * sp.sin(x)
c3x = sp.cos(x) * sp.cos(x) * sp.cos(x)
f   = sp.exp(x)/sp.sqrt(s3x + c3x)
df  = sp.diff(f, x)

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Function, f: " + bcolors.ENDC)
print("")
sp.pprint(f)
print("")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Derivative, f': " + bcolors.ENDC)
print("")
sp.pprint(df)
print("")

#EOF