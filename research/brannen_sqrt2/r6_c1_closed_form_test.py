"""
R6.8: Czy c1 = 0.72538 (asymetria deficit/excess) ma prosta postac zamknieta?

c1 wynika z eta(delta) = A_tail(1+delta)/|delta|:
  eta_exc(delta) - eta_def(delta) = c1 * delta  (do rzedu linear)

Hipoteza: c1 to klasyczny ulamek zwiazany z K=2/3 lub Z3 symetria.
"""
import math

c1 = 0.72538

candidates = {
    "ln(3)/4":              math.log(3)/4,
    "1 - ln(3)/4":          1 - math.log(3)/4,
    "1 - ln(3)/4 + e-5":    1 - math.log(3)/4 - 1.3e-4,
    "3/4 - 1/40":           3/4 - 1/40,
    "3/4 * (1-1/32)":       3/4 * (1 - 1/32),
    "2/3 + 1/17":           2/3 + 1/17,
    "(pi - 1)/e":           (math.pi-1)/math.e,
    "sqrt(phi) - 1/2":      math.sqrt(1.6180339887) - 0.5,
    "(1 + 1/phi)/2":        (1 + 1/1.6180339887)/2,
    "phi/sqrt(5)":          1.6180339887/math.sqrt(5),
    "2/e":                  2/math.e,
    "1 - 1/e^(3/2)":        1 - math.exp(-1.5),
    "(phi-1)/phi * 2":      (1.6180339887-1)/1.6180339887 * 2,
    "pi / e / sqrt(2)":     math.pi/math.e/math.sqrt(2),
    "sqrt(2) - 1/sqrt(2)":  math.sqrt(2) - 1/math.sqrt(2),
    "arctan(pi/2)":         math.atan(math.pi/2),
    "1 - 1/pi":             1 - 1/math.pi,
    "pi^2/8 - 1/pi":        math.pi**2/8 - 1/math.pi,
    "3 ln(2) / pi":         3*math.log(2)/math.pi,
    "(e-1)/sqrt(e)":        (math.e-1)/math.sqrt(math.e),
    "log(3)/log(2)/2":      math.log(3)/math.log(2)/2,
    "3/4 - 1/sqrt(3)/100":  0.75 - 1/math.sqrt(3)/100,
    "cos(pi/5)/sqrt(2)":    math.cos(math.pi/5)/math.sqrt(2),
    "1/sqrt(2) + 1/50":     1/math.sqrt(2) + 1/50,
    "gamma(1/3)/gamma(2/3)/2": 2.6789385347/1.3541179394/2,  # incomplete - rough
    "4*(3/4)^3":            4*(0.75)**3,
    "7/(3*pi)":             7/(3*math.pi),
    "sqrt(3)/(pi-phi)":     math.sqrt(3)/(math.pi-1.618),
    "2/3 * (1 + 1/(4*sqrt(pi)))": (2/3)*(1 + 1/(4*math.sqrt(math.pi))),
    "2/3 + 1/(16 ln 3)":    2/3 + 1/(16*math.log(3)),
    "(K^-1) * (ln2)":       1.5 * math.log(2),  # 3/2 * ln(2) = 1.04
    "1/ln(e*phi)":          1/math.log(math.e*1.618),
}

print("c1 numerical = {:.8f}".format(c1))
print()
header = "{:30s} {:>15s} {:>14s} {:>10s}".format("Candidate", "Value", "Diff", "Rel err")
print(header)
print("-" * len(header))
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-c1)):
    diff = val - c1
    rel = abs(diff)/c1
    print("{:30s} {:15.8f} {:+14.4e} {:10.2e}".format(name, val, diff, rel))

# Szczegolowy test: c1 vs 1 - ln(3)/4
print()
print("=== Szczegolowa analiza: 1 - ln(3)/4 hipoteza ===")
target = 1 - math.log(3)/4
print("  c1                = {:.10f}".format(c1))
print("  1 - ln(3)/4       = {:.10f}".format(target))
print("  c1 - (1-ln3/4)    = {:+.6e}".format(c1 - target))
print("  |diff|/c1         = {:.2e}".format(abs(c1-target)/c1))
print()
print("  4*(1-c1)          = {:.6f}".format(4*(1-c1)))
print("  ln(3)             = {:.6f}".format(math.log(3)))
print("  4(1-c1) - ln(3)   = {:+.6e}".format(4*(1-c1) - math.log(3)))

# Polepszenie: c1 pochodzi z calki/asymptotyki ODE
# Moze bardziej naturalna forma: pokazac ze 4(1-c1) = ln(3) exactly
# lub: exp(4(1-c1)) = 3 exactly
print()
print("  exp(4*(1-c1))     = {:.6f}   (hipoteza: = 3)".format(math.exp(4*(1-c1))))
print("  exp(ln 3)         = {:.6f}".format(math.exp(math.log(3))))
print("  diff              = {:+.6e}".format(math.exp(4*(1-c1)) - 3))
