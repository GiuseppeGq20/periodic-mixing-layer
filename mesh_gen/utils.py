import gmsh
import math

def calc_progressin(L,l0,N):
    """Evaluate the parameter a in geometric series such that:

    L = l_{0}\\sum_{i=0}^{N} a^{i}
    
    It uses a secant method to  find the root of the function.

    Args:
        L - Total length
        l0 - minimum length
        N - number of subdivision

    Returns:
        a - geometric series growth factor
    """
    n = N -1

    def f(p):
        #f = (p**(N-1))*(1+p) - (L/l0)
        f = (L/l0) - (1- p**(N+1))/(1-p)
        return f

    x1 = 1.1
    x2 = 1.3
    tol = 1e-4
    err = 1e2
    iter_max = 1000
    i = 0
    while (err > tol and i < iter_max):
        f1 = f(x1)
        f2 = f(x2)
        x_new = x2 - f2*(x2-x1)/(f2-f1)
        x1 = x2
        x2 = x_new
        err = abs(f(x_new))
        i +=1
    
    print("iter: ", i)
    return x2

def vertical_points(h: list[float], x_coord: float,lc=0)->list[int]:
    """generate a list of gmsh point located on the line z=0, x = x_coord with the gmsh geometric kernel

    Args:
        h (list): y coordinate of the points
        x_coord (float): x coordinate of the points
        lc (int, optional): gmsh mesh element sizing. Defaults to 0.

    Returns:
        p: list containing the generated point tags
    """
    p = []
    for i in range(len(h)):
        p.append(
            gmsh.model.geo.addPoint(x_coord, h[i], 0, lc)
        )
    return p

def horizontal_lines(p1: list[int],p2: list[int]) -> list[int]:
    """return a list of the gmsh lines generated with the point list provided using the 
    gmsh geometric kernel

    Args:
        p1 (list[int]): list with the tags of the first point of the line
        p2 (list[int]): list with the tags of the second point of the line

    Returns:
        l (list[int]): list containing the tags of the generated lines
    """
    if len(p1) != len(p2): raise(ValueError)
    l = []
    for i in range(len(p1)):
        l.append(gmsh.model.geo.addLine(p1[i],p2[i]))
    return l

def vertical_lines(p1: list[int])->list[int]:
    """returns a list of the gmsh line generated with the point tag list provided, 
    it uses the gmsh geometric kernel.
    p1 points must lie on the same line

    Args:
        p1 (list[int]): list of input point tags

    Returns:
        l (list[int]): list containing the tags of the generated lines
    """
    if (len(p1) <= 1):
        raise(ValueError)
    l  = []

    for i in range(1,len(p1)):
        l.append(gmsh.model.geo.addLine(p1[i],p1[i-1]))
    return l

def quad_loops(
        v_left_lines: list[int],
        v_right_lines: list[int],
        bottom_lines: list[int],
        top_lines: list[int]) -> list[int]:
    """Returns a list of the curve loop tags of the curve loop generated.
    It generates quadrilateral following a counterclowise path from left, bottom, right and top side.
    Lines list provided must be in such orded. It uses the gmsh geometric kernel

    Args:
        v_left_lines (list[int]): list of left side tags
        v_right_lines (list[int]): list of right side tags
        bottom_lines (list[int]): list of bottom side tags

    Returns:
        q_loops list[int]: list of quad loop tags
    """
    vl,vr = v_left_lines, v_right_lines
    bl,tl = bottom_lines, top_lines
    if len(vl) != len(vl) and len(vl) != len(bl) and len(bl) != len(tl): raise(ValueError)

    N = len(bl)
    c_loop = []
    for i in range(N):
        c_loop.append(
            gmsh.model.geo.addCurveLoop([vl[i], bl[i], vr[i], tl[i]], -1,reorient=True)
        )
    return c_loop
