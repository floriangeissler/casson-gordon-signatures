"""
Copyright (c) 2025 Florian Geissler

This software is released under the MIT License.
See the LICENSE file for details.

A script to compute Casson-Gordon sigma signatures of two-bridge knots
using the formula from [CG1986, p. 188]. The formula expresses the
signatures in terms of a weighted count of lattice points in a planar
triangle determined by the knot.

References:
    [CG1986] Casson, A.J. and Gordon, C. McA., "Cobordism of Classical Knots",
    Progress in Mathematics, Vol. 62, pp. 181–199, 1986.
"""

import math
from fractions import Fraction
from typing import Dict

try:
    from matplotlib.ticker import MaxNLocator
    from matplotlib import pyplot
    _is_matplotlib_available = True
except ImportError:
    _is_matplotlib_available = False
    print("Warning: Matplotlib not found. Plotting will be disabled.")


def ensure_positive_integer(value: int) -> int:
    """
    Ensures input is an integer greater than zero.
    """
    if not isinstance(value, int) or value <= 0:
        raise ValueError(f"Input must be a positive int value, but received: {value} of type {type(value)}.")
    return value


def ensure_positive_fraction(value: int | float | Fraction) -> Fraction:
    """
    Ensures input greater than zero and returns it as a `Fraction` data type.
    """
    if not isinstance(value, (int, float, Fraction)) or not value > 0:
        raise ValueError(f"Input must be a positive int, float, or Fraction, but got {value} of type {type(value)}.")
    return Fraction(value)



def compute_weighted_vertex_number(delta_x: Fraction, delta_y: Fraction, debug: bool = False, plot: bool = False) -> Fraction:
    """
    Computes the weighted vertex number for a triangle ∆(x,y) following Casson-Gordon [CG1986, p. 188].
    The parameters x, y are determined by the Schubert normal form of 2-bridge knot and parameters r and m, 
    following [CG1986, p. 187]. The sum of integer lattice points is calculated with the following weights:
        - interior points: 1
        - boundary points: 1/2
        - vertices B and C: 1/4
        - vertex A (0,0): 0 (not counted)

    The triangle ∆(x,y) is depicted below:

              (x,y)
                C
               /|
              / |   
             /  |     
            /   |        
           /    |
          /     |
         /      |
        /_______| 
       A        B
     (0,0)    (x,0)


    Args:
        delta_x (Fraction): The x-coordinate of ∆(x,y).
        delta_y (Fraction): The y-coordinate of ∆(x,y) .
        debug (bool): If True, prints debug information.
        plot (bool): If True, generates a plot of the triangle and lattice points.

    Returns:
        Fraction: The weighted vertex number.

    Raises:
        ValueError: If delta_x or delta_y are not > 0.

    References:
        [CG1986] Casson, A.J. and Gordon, C. McA., "Cobordism of Classical Knots",
        Progress in Mathematics, Vol. 62, pp. 181–199, 1986.
    """

    delta_x = ensure_positive_fraction(delta_x)
    delta_y = ensure_positive_fraction(delta_y)
    delta_zero = Fraction(0)

    # Define lattice point weights as per [CG1986, p. 188]
    interior_point_weight = Fraction(1)
    boundary_point_weight = Fraction(1, 2)
    vertex_point_weight = Fraction(1, 4)
    
    # Initialize lists for lattice point counting
    interior_points = []
    boundary_points = []
    vertex_points = []
 

    
    # Define the triangle's sloped edge equation
    slope = delta_y / delta_x

    # The nested loops iterate through all integer lattice points (px, py) which lie in ∆(x,y).
    # The `ceil` function handles fractional inputs. We add 1 to the ceiling to ensure 
    # vertices B and C are included if their coordinates are integers.
    for px in range(math.ceil(delta_x) + 1):
        for py in range(math.ceil(delta_y) + 1):
            
            # Vertex A (0,0) is excluded from the count.
            if px == delta_zero and py == delta_zero:
                # The point is on the sloped edge. We break the inner loop,
                # since all subsequent py values for this px will be above.
                break

            # Integer lattice points (px, py) lie inside or on a boundary of the triangle if
            # ( px <= delta_x ) and ( py <= slope * px )
        
            if px > delta_x:
                # The point is outside the vertical edge. We break the inner loop,
                # since all subsequent py values for this px will also be outside.
                # We need this check in cases when delta_x is not an integer.
                break

            if py > px * slope:
                # The point is above the sloped edge. We break the inner loop,
                # since all subsequent py values for this px will also be above.
                break

            # Classify boundary points and vertex points
            is_on_boundary_ab = (py == delta_zero) and (px <= delta_x)
            is_on_boundary_bc = (px == delta_x) and (py <= delta_y)
            is_on_boundary_ac = (py == px * slope) and (py > delta_zero)
            is_vertex_b = (px == delta_x) and (py == delta_zero)
            is_vertex_c = (px == delta_x) and (py == delta_y)
            
            # Add points to their respective lists
            if is_vertex_b or is_vertex_c:
                # Point is vertex B or C
                vertex_points.append((px, py))
            elif is_on_boundary_bc or is_on_boundary_ab or is_on_boundary_ac:
                # Point is boundary point other than B or C
                boundary_points.append((px, py))
            else:
                # All other points are interior points 
                interior_points.append((px, py))

    
    # Compute the weighted vertex number `int ∆(x, y)`
    weighted_vertex_number = (
          (len(interior_points) * interior_point_weight)
        + (len(boundary_points) * boundary_point_weight)
        + (len(vertex_points) * vertex_point_weight)
    )

    if debug:
        print(f"∆({delta_x}, {delta_y}) vertices are A({delta_zero},{delta_zero}), B({delta_x},{delta_zero}), C({delta_x},{delta_y}) with "
              f"slope = ({slope}) ≈ {round(float(slope), 3)}")
        print(f'Interior Points:  count = {len(interior_points)}\tList:', interior_points)
        print(f'Boundary Points:  count = {len(boundary_points)}\tList:', boundary_points)
        print(f'Vertex Points:    count = {len(vertex_points)}  \tList:', vertex_points)

        area = (delta_x * delta_y / 2)
        print(f"Weighted Vertex Number = {weighted_vertex_number}",
              f"≈ {round(float(weighted_vertex_number), 3)}" if weighted_vertex_number.denominator != 1 else "",
              "is equal to" if (area == weighted_vertex_number) else "is NOT EQUAL to",
              f"area={area}",
              f"(≈ {round(float(area), 3)})" if area.denominator != 1 else ""
              )

    # Plot the diagram for the triangle ∆(x, y) (optional)
    if plot:
        plot_diagram(interior_points, boundary_points, vertex_points, delta_x, delta_y)

    return weighted_vertex_number


def compute_cgsigma(p: int, q: int, m: int, debug: bool = False, plot: bool = False) -> Dict[int, Fraction]:
    """
    Computes the Casson-Gordon sigma value for parameters p, q, m and each r in {1, ..., m-1}
    using the formula from [CG1986 p.188]: `cgsigma = 4 * ( area ∆(x,y) - int ∆(x,y) )`

    Args:
        p (int): The parameter p.
        q (int): The parameter q.
        m (int): The parameter m.
        debug (bool): If True, prints debug information.
        plot (bool): If True, generates a plot of the triangles for each r.

    Returns:
        dict: A dictionary of {r: cgsigma} for all r in {1, ..., m-1}.
    """
    p = ensure_positive_integer(p)
    q = ensure_positive_integer(q)
    m = ensure_positive_integer(m)

    cgsig_values = {}
    for r in range(1, m):

        delta_x = Fraction(p * r, m)
        delta_y = Fraction(q * r, m)
        
        # area ∆(x, y) = 1/2 * base * height = (delta_x * delta_y)/2
        area = Fraction(delta_x * delta_y, 2)

        # int ∆(x, y) is the weighted vertex number
        weighted_vertex_number = compute_weighted_vertex_number(delta_x, delta_y, debug=debug, plot=plot)

        # Casson-Gordon signature formula: 4 * (area ∆(x, y) - int ∆(x, y))
        cgsig_value = 4 * (area - weighted_vertex_number)
        cgsig_values[r] = int(cgsig_value) if cgsig_value.denominator == 1 else cgsig_value

        if debug:
            print(f"p={p}, q={q}, m={m}, r={r}, x={delta_x}, y={delta_y}")
            print(f"WVN({r})={weighted_vertex_number}, area({r})={area}")
            print(f"cgsigma({r})={cgsig_value}")
            print("---------------------")
            
    return cgsig_values


def plot_diagram(interior_points, boundary_points, vertex_points, delta_x, delta_y):
    """
    Visualizes the triangle ∆(x,y) and lattice points using matplotlib.

    Parameters:
        interior_points (list of tuples): Coordinates of interior lattice points.
        boundary_points (list of tuples): Coordinates of boundary lattice points.
        vertex_points (list of tuples): Coordinates of vertex points.
        delta_x, delta_y (Fraction): x- and y-coordinates for the triangle vertices.
    """

    if not _is_matplotlib_available:
        print("Plotting skipped: Matplotlib not installed.")
        return

    # Skip plotting unreadable diagrams.
    if len(interior_points) > 200:
        print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
              "Skipping plot for ∆(x,y) with over 200 interior_points."
              "\n- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
        return

    # -- Plotting the triangle
    fig, ax = pyplot.subplots()

    # -- Initialize variables for vertex coordinates
    Ax, Ay = int(0), int(0)
    Bx, By = float(delta_x), int(0)
    Cx, Cy = float(delta_x), float(delta_y)

    # -- Plot points
    # See https://matplotlib.org/stable/gallery/color/named_colors.html for different color options.
    for point in interior_points:
        ax.plot(point[0], point[1],  marker='o', color='darkorange', alpha=0.5)
    for point in boundary_points:
        ax.plot(point[0], point[1], marker='D', color='green', alpha=0.5)
    for point in vertex_points:
        ax.plot(point[0], point[1],   marker='s', color='royalblue', alpha=0.75)


    # -- Plot line segments [A,B], [A,C], [B,C]
    ax.plot([Ax, Bx, Cx, Ax], [Ay, By, Cy, Ay], color='black', alpha=0.8, linewidth=2)

    # -- Label the vertex points with an offset
    label_offset = max(Bx, Cy) * 0.02
    ax.text(Ax - label_offset, Ay - label_offset, 'A', fontsize=12, ha='right', va='top', fontweight='bold')
    ax.text(Bx + label_offset, By - label_offset, 'B', fontsize=12, ha='left', va='top', fontweight='bold')
    ax.text(Cx + label_offset, Cy + label_offset, 'C', fontsize=12, ha='left', va='bottom', fontweight='bold')

    # -- Set axis limits for better visibility
    margin_x = max(1, Bx * 0.2)
    margin_y = max(1, Cy * 0.2)
    ax.set_xlim([-margin_x, math.ceil(Cx) + margin_x])
    ax.set_ylim([-margin_y, math.ceil(Cy) + margin_y])
    
    # -- Set the x and y axis ticks at integer values
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    # -- Title and labels
    ax.set_title(f'Triangle ∆({delta_x}, {delta_y}) and Lattice Points')
    ax.set_xlabel('x-axis')
    ax.set_ylabel('y-axis')
    
    # -- Show grid and plot
    ax.grid(True)
    pyplot.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Compute CGSigma values for given p, q, m."
    )
    parser.add_argument("p", type=int, help="The value for p (must be a positive integer)")
    parser.add_argument("q", type=int, help="The value for q (must be a positive integer)")
    parser.add_argument("m", type=int, help="The value for m (must be a positive integer)")
    parser.add_argument("--debug", action="store_true", help="Enable debug print statements.")
    parser.add_argument("--plot", action="store_true", help="Generate plots for each triangle.")
    
    args = parser.parse_args()

    try:
        p_val = ensure_positive_integer(args.p)
        q_val = ensure_positive_integer(args.q)
        m_val = ensure_positive_integer(args.m)
        
        results = compute_cgsigma(p_val, q_val, m_val, args.debug, args.plot)
        print("\n--- CGSigma Values ---")
        for r, value in results.items():
            print(f"r={r}: {value}")

    except ValueError as e:
        print(f"Error: {e}")
        parser.print_help()