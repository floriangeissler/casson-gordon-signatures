# CGSigma: Casson-Gordon Signature Calculator

A Python script to compute Casson-Gordon sigma signatures of two-bridge knots using the formula from [CG1986, p. 188]. The formula expresses the signatures in terms of a weighted count of lattice points in a planar triangle determined by the knot.

The script can optionally generate plots of these triangles and their corresponding lattice points using matplotlib.



## Usage: Command Line 
```bash
python cgsigma.py <p> <q> <m> [--debug] [--plot]
```

#### Example: Compute CGSigma values
```bash
python cgsigma.py 9 2 3
```

#### Example: Compute CGSigma values with debug information and plots
```bash
python cgsigma.py 9 2 3 --debug --plot
```

## Usage: Jupyter Notebook

#### Example: Compute CGSigma values
```python
from cgsigma import compute_cgsigma

p,q,m = 9,2,3

cgsigma_values = compute_cgsigma(int(p), int(q), int(m), debug=False, plot=False)

print(f"CGSigma values for all r:", cgsigma_values)

r=2
print(f"CGSigma values for r={r}:", cgsigma_values[r])
```

#### Example: Compute Weighted Vertex Number with debug information and plots
```python
from cgsigma import compute_weighted_vertex_number
from fractions import Fraction

delta_x=Fraction(int(9),int(1))
delta_y=int(6)

wvn = compute_weighted_vertex_number(delta_x, delta_y, debug=True, plot=True)

print(f"Weighted Vertex Number:", wvn)
```

## References
[CG1986] Casson, A.J. and Gordon, C. McA., "Cobordism of Classical Knots", Progress in Mathematics, Vol. 62, pp. 181–199, 1986.

## License
This code is licensed under the MIT License. See the LICENSE file for details.