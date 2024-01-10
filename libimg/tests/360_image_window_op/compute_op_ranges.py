#!/usr/bin/python3
# Last edited on 2023-09-10 19:43:11 by stolfi

# Computes ranges and other information about the "smoothed" 3x4 window operators.

from math import sqrt, hypot, inf, floor;
from sys import exit, stderr as err;
import re

n = 9  # Number of basis elements.
ones = [[1,1,1],[1,1,1],[1,1,1]]

def test_basis(smooth):
  err.write(f"### test_basis, smooth = {smooth} ##################################################\n")
  B, W = build_basis(smooth)
  
  test_basis_single(B, W, "ones", [[+1,+1,+1],[+1,+1,+1],[+1,+1,+1]])
  test_basis_single(B, W, "X",    [[-1,00,+1],[-1,00,+1],[-1,00,+1]])
  test_basis_single(B, W, "Y",    [[-1,-1,-1],[00,00,00],[+1,+1,+1]])
  test_basis_single(B, W, "X^2",  [[+1,00,+1],[+1,00,+1],[+1,00,+1]])
  test_basis_single(B, W, "Y^2",  [[+1,+1,+1],[00,00,00],[+1,+1,+1]])
  test_basis_single(B, W, "XY",   [[+1,00,-1],[00,00,00],[-1,00,+1]])
  test_basis_single(B, W, "1 + X^2 + Y^2", [[+3,+2,+3],[+2,+1,+2],[+3,+2,+3]])
  # ----------------------------------------------------------------------

def test_basis_single(B, W, name, Z):
  C = decomp_window(Z, B, W)
  err.write("%-12s %a  ->  %a\n" % (name + ":", Z,C))
  # ----------------------------------------------------------------------

def test_ops():
  err.write(f"### test_ops ##################################################\n")
  Bu, Wu = build_unsmooth_basis()
  Bs, Ws = build_smooth_basis()
  
  # Samples from a constant function:
  Zname = "constant {Z = 7}"
  Z = [[7,7,7],[7,7,7],[7,7,7]]
  OPue = {
    "f": 7, "fx": 0, "fy": 0, "fxx": 0, "fxy": 0, "fyy": 0,
    "fxxy": 0, "fxyy": 0, "fxxyy": 0,
    "laplacian": 0, "orthicity": 0, "elongation": 0, "saddleness": 0,
    "const1": 7, "constdev": 0,
    "lin1": 7, "linX": 0, "linY": 0, "lindev": 0,
    "quad1": 7, "quadX": 0, "quadY": 0, "quadX2": 0, "quadXY": 0, "quadY2": 0, "quaddev": 0,
    "average": 7, "deviation": 0,
  }
  OPse = OPue
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)
  
  # Samples from a linear function:
  Zname = "linear {Z = 2 + 3*X + 4*Y}"
  Z = [[-5,-2,+1],[-1,+2,+5],[+3,+6,+9]]
  OPue = {
    "f": 2, "fx": 3, "fy": 4, "fxx": 0, "fxy": 0, "fyy": 0,
    "fxxy": 0, "fxyy": 0, "fxxyy": 0,
    "laplacian": 0, "orthicity": 0, "elongation": 0, "saddleness": 0,
    "const1": 2, "constdev": 5*sqrt(0.5),
    "lin1": 2, "linX": 3, "linY": 4, "lindev": 0,
    "quad1": 2, "quadX": 3, "quadY": 4, "quadX2": 0, "quadXY": 0, "quadY2": 0, "quaddev": 0,
    "average": 2, "deviation": 5*sqrt(0.5),
  }
  OPse = OPue.copy()
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)
  
  # Samples from a quadratic function:
  Zname = "quadratic {Z = 2 - 3*X + 4*Y - 5*X^2 + 6*X*Y - 7*Y^2}"
  Z = [[-5,-9,-23],[0,+2,-6],[-9,-1,-3]]
  OPue = {
    "f": 2, "fx": -3, "fy": 4, "fxx": -10, "fxy": 6, "fyy": -14,
    "fxxy": 0, "fxyy": 0, "fxxyy": 0,
    "laplacian": -24, "orthicity": +4, "elongation": 4*sqrt(10), "saddleness": 244,
    "const1": 2, "constdev": sqrt(76),
    "lin1": 2, "linX": -3, "linY": +4, "lindev": sqrt(63.5),
    "quad1": 2, "quadX": -3, "quadY": 4, "quadX2": -5, "quadXY": +6, "quadY2": -7, "quaddev": 0,
    "average": -4, "deviation": 2*sqrt(10),
  }
  OPse = OPue.copy()
  OPse["const1"] = -4
  OPse["constdev"] = 2*sqrt(10)
  OPse["lin1"] = -4
  OPse["lindev"] = sqrt(27.5)
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)
  
  # Samples from a cubic type 1 function:
  Zname = "cubic type 1 {Z = 3*X^2*Y + 5*X*Y^2}"
  # Should be recognized by unsmoothed ops, not by smoothed ones.
  Z = [[-8,0,+2],[0,0,0],[-2,0,+8]]
  OPue = {
    "f": 0, "fx": 0, "fy": 0, "fxx": 0, "fxy": 0, "fyy": 0,
    "fxxy": 6, "fxyy": 10, "fxxyy": 0,
    "laplacian": 0, "orthicity": 0, "elongation": 0, "saddleness": 0,
    "const1": 0, "constdev": sqrt(17/2),
    "lin1": 0, "linX": 0, "linY": 0, "lindev": sqrt(17/2),
    "quad1": 0, "quadX": 0, "quadY": 0, "quadX2": 0, "quadXY": 0, "quadY2": 0, "quaddev": sqrt(17/2),
    "average": 0, "deviation": sqrt(17/2), 
  }
  OPse = OPue.copy()
  OPse["fx"] = 2.5
  OPse["fy"] = 1.5
  OPse["linX"] = 2.5
  OPse["linY"] = 1.5
  OPse["lindev"] = 2
  OPse["quadX"] = 2.5
  OPse["quadY"] = 1.5
  OPse["quaddev"] = 2
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)

  # Samples from the cubic  type 2:
  Zname = "cubic type 2 {Z = 3*X*(X^2 - 2*Y^2) + 5*Y*(Y^2 - 2*X^2)}"
  # Should be recognized by smoothed ops, not by unsmoothed ones.
  Z = [[+8,-5,+2],[-3,0,+3],[-2,+5,-8]]
  OPue = {
    "f": 0, "fx": 3, "fy": 5, "fxx": 0, "fxy": 0, "fyy": 0,
    "fxxy": -20, "fxyy": -12, "fxxyy": 0,
    "laplacian": 0, "orthicity": 0, "elongation": 0, "saddleness": 0,
    "const1": 0, "constdev": sqrt(17),
    "lin1": 0, "linX": 3, "linY": 5, "lindev": sqrt(34),
    "quad1": 0, "quadX": 3, "quadY": 5, "quadX2": 0, "quadXY": 0, "quadY2": 0, "quaddev": sqrt(34),
    "average": 0, "deviation": sqrt(17), 
  }
  OPse = OPue.copy()
  OPse["fx"] = 0
  OPse["fy"] = 0
  OPse["linX"] = 0
  OPse["linY"] = 0
  OPse["lindev"] = sqrt(17)
  OPse["quadX"] = 0
  OPse["quadY"] = 0
  OPse["quaddev"] = sqrt(17)
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)
  
  # Samples from a special quartic:
  Zname = "quartic {Z = 2 + 3*Y - 5*X^2 + 6*X^2*Y^2}"
  # Should be recognized by smoothed ops?
  Z = [[0,-1,0],[-3,2,-3],[6,5,6]]
  OPue = {
    "f": 2, "fx": 0, "fy": 3, "fxx": -10, "fxy": 0, "fyy": 0,
    "fxxy": 0, "fxyy": 0, "fxxyy": 24,
    "laplacian": -10, "orthicity": -10, "elongation": 10, "saddleness": 0,
    "const1": 2, "constdev": sqrt(11),
    "lin1": 2, "linX": 0, "linY": 3, "lindev": sqrt(8),
    "quad1": 2, "quadX": 0, "quadY": 3, "quadX2": -5, "quadXY": 0, "quadY2": 0, "quaddev": sqrt(12.5),
    "average": 1, "deviation": sqrt(229/9), 
  }
  OPse = OPue.copy()
  OPse["f"] = 1
  OPse["fx"] = None
  OPse["fy"] = None
  OPse["fxx"] = None
  OPse["fyy"] = None
  OPse["const1"] = 1
  OPse["constdev"] = sqrt(10)
  OPse["lin1"] = 1
  OPse["linX"] = None
  OPse["linY"] = None
  OPse["lindev"] = None
  OPse["quad1"] = 1
  OPse["quadX"] = None
  OPse["quadY"] = None
  OPse["quadX2"] = None
  OPse["quadY2"] = None
  OPse["quaddev"] = None
  
  test_ops_one(Zname, Z, Bu, Wu, OPue, Bs, Ws, OPse)
  # ----------------------------------------------------------------------
  
def test_ops_one(name, Z, Bu, Wu, OPue, Bs, Ws, OPse):
  # Tests the unsmoothed and smoothed operators on {Z}, given their bases
  # and the expected results.
  err.write("-"*70 + "\n")
  err.write(f"testing {name}\n")
  err.write("%-10s %21s %2s   %21s\n" % ("op","unsmooth","","smooth"))
  err.write("%-10s %21s %2s   %21s\n" % ("-"*10,"-"*21,"","-"*21))
  OPu = compute_unsmooth_ops(Z, Bu, Wu)
  OPua = compute_avg_dev_ops(Z, Ws, OPu) # YES, {Ws} not {Wu} here.
  OPu |= OPua; 
  OPs = compute_smooth_ops(Z, Bs, Ws);
  OPsa = compute_avg_dev_ops(Z, Ws, OPs)
  OPs |= OPsa
  print_ops4(OPu,OPue,OPs,OPse)
  err.write("-"*70 + "\n")
  # ----------------------------------------------------------------------

def test_enum_windows(smooth):

  err.write(f"### test_enum_windows, smooth = {smooth} ##################################################\n")
  B, W = build_basis(smooth)

  # Pixel values:
  Z = [[None,None,None],[None,None,None],[None,None,None]]

  Min = {}; ZMin = {}
  Max = {}; ZMax = {}
  
  q = 3 # Levels per element.

  # Enumerate all windows with {q} values per element:
  for b in range(q**n):

    # Convert {b} to a 0-1 window:
    bb = b;
    for Y in -1,0,+1:
      for X in -1,0,+1:
        Z[Y+1][X+1] = (bb % q)/(q - 1)
        bb = bb // q
    assert bb == 0
    OP = compute_ops(Z, smooth, B, W)
    for key in OP.keys():
      accum(Min, ZMin, Max, ZMax, key, OP[key], Z)
  print_ranges(Min, ZMin, Max, ZMax)
  # ----------------------------------------------------------------------

def build_basis(smooth):
  if smooth:
    B, W = build_smooth_basis()
  else:
    B, W = build_unsmooth_basis()
  return B, W
  # ----------------------------------------------------------------------

def compute_ops(Z, smooth, B, W):
  # {B,W} must be appropriate for {smooth}.
  if smooth:
    OP = compute_smooth_ops(Z, B, W)
  else:
    OP = compute_unsmooth_ops(Z, B, W)
  return OP
  # ----------------------------------------------------------------------

def unsmooth_weight_window():
  # The weight window to be used for dot products in unsmooth operators etc.
  # The weights are all 1.0
   
  W = [[0,0,0],[0,0,0],[0,0,0]]
  for i in range(3):
    for j in range(3):
      W[i][j] = 1.0
  assert abs(dot(ones, ones, W) - 9) < 1.0e-13;
  return W
  # ----------------------------------------------------------------------

def build_unsmooth_basis():
  # Returns the basis {B} of the unsmooth operators, as a list of {n} 3x3
  # arrays, and the corresponding weight window {W}. Also checks that
  # they are orthogonal.

  W = unsmooth_weight_window()

  assert n == 9
  # In the following comments, {X} and {Y} are the windows with 
  # {Xij = j-1} and {Yij = i-1}.
  B = ( \
    ((00,00,00),(00,+1,00),(00,00,00)),  # {(1-X^2)*(1-Y^2)}.
    ((00,00,00),(-1,00,+1),(00,00,00)),  # {X*(1-Y^2)}.
    ((00,-1,00),(00,00,00),(00,+1,00)),  # {Y*(1-X^2)}.
    ((00,00,00),(+1,00,+1),(00,00,00)),  # {X^2*(1 - Y^2)}.
    ((+1,00,-1),(00,00,00),(-1,00,+1)),  # {X*Y}.
    ((00,+1,00),(00,00,00),(00,+1,00)),  # {Y^2*(1 - X^2)}.
    
    ((-1,00,-1),(00,00,00),(+1,00,+1)),  # {X^2*Y}.
    ((-1,00,+1),(00,00,00),(-1,00,+1)),  # {X*Y^2}.
    ((+1,00,+1),(00,00,00),(+1,00,+1)),  # {X^2*Y^2}.
  )
  #
  # Note that
  #
  #  {ones = B[0]+B[3]+B[5]+B[8]}
  #  {X = B[1]+B[6]}
  #  {Y = B[2]+B[7]}
  #  {X^2 = B[3]+B[8]}
  #  {Y^2 = B[5]+B[8}}
  #
  
  # Check orthogonality:
  for r in range(n):
    for s in range(n):
      Mrs = dot(B[r], B[s], W)
      if r == s:
        if Mrs <= 0:
          err.write(f"** dot <B{r}|B{r}> = {Mrs}\n"); assert False
      else:
        if Mrs != 0:
          err.write(f"** dot <B{r}|B{s}> = {Mrs}\n"); assert False
  return B, W
  # ----------------------------------------------------------------------

def compute_unsmooth_ops(Z, B, W):
  # Computes all the 3x3 unsmooth operators from the window {Z[0..2][0..2]}.
  # Returns them as a dict.

  C = decomp_window(Z, B, W)
  OP = {}
  
  # Differential operators: computes the estimates of the differential
  # operators at the center of the window, by the simplest unbiased
  # discrete difference formula.
  #
  # This is equivalent to assuming that the "cardinal" elements of the
  # weight window {W01,W21,W10,W12} are implicitly multiplied by an
  # infinitesimal {eps}, while the "corner" elements {W00,W02,W20,W22}
  # are implicitly multiplied by {eps^2}. Thus the the coefficients
  # {C[1],C[2],C[3],C[5]} of the "cardinal" basis elements are
  # implicitly multiplied by {eps} while the coefficients
  # {C[4],C[6],C[7],C[8]} are implicitly multiplied by {eps^2}

  f = C[0];                              OP["f"] = f
  fx = C[1];                             OP["fx"] = fx
  fy = C[2];                             OP["fy"] = fy
  fxx = 2*(C[3] - C[0]);                 OP["fxx"] = fxx
  fxy = C[4];                            OP["fxy"] = fxy
  fyy = 2*(C[5] - C[0]);                 OP["fyy"] = fyy
  fxyy = 2*C[6] - 2*C[1];                OP["fxyy"] = fxyy
  fxxy = 2*C[7] - 2*C[2];                OP["fxxy"] = fxxy
  fxxyy = 4*(C[0]-C[3]-C[5]+C[8]);       OP["fxxyy"] = fxxyy

  laplacian = fxx+fyy;                   OP["laplacian"] = laplacian
  orthicity = fxx-fyy;                   OP["orthicity"] = orthicity
  elongation = hypot(2*fxy, orthicity);  OP["elongation"] = elongation
  saddleness = 2*fxx*fyy - fxy*fxy;      OP["saddleness"] = saddleness

  # Computes the zero-order Taylor approximation {const1}
  # to {Z[0..2][0..2]} based on the unsmooth differential operators.
  # Namely, the coefficient {const1} is just {f} namely the central pixel:
  const1 = C[0];       OP["const1"] = const1
  
  # Computes the Taylor approximation {lin1 + linX*X + linY*Y} 
  # to {Z[0..2][0..2]} based on the unsmooth differential operators.
  # Namely, the coefficient {lin1} is just the {f}, and the coefficients
  # {linX,linY} are just {fx} and {fy}, respectively. 
  
  lin1 = f;  OP["lin1"] = lin1
  linX = fx;     OP["linX"] = linX
  linY = fy;     OP["linY"] = linY

  # Computes the quadratic Taylor approximation to {Z[0..2][0..2]}
  # {quad1 + quadX*X + quadY*Y + quadX2*X^2 + quadXY*X*Y + quadY2*Y^2}
  # based on the unsmooth differential operators. Namely,
  # {quad1,quadX,quadY,quadXY} are {f,fx,fy,fxy}, while
  # {quadX2} and {quadY2} are {fxx/2} and {fyy/2}.
  
  quad1 = f;        OP["quad1"] = quad1   # Same as the "f" oeprator
  quadX = fx;       OP["quadX"] = quadX   # Same as "fx".
  quadY = fy;       OP["quadY"] = quadY   # Same as "fy".
  quadX2 = fxx/2;   OP["quadX2"] = quadX2 # same as "fxx"/2.
  quadXY = fxy;     OP["quadXY"] = quadXY # same as "fxy".
  quadY2 = fyy/2;   OP["quadY2"] = quadY2 # same as "fyy"/2.

  return OP
  # ----------------------------------------------------------------------

def smooth_weight_window():
  # The weight window to be used for dot products in smooth operators. etc.
  # The weights are normalized to add to 1.
  W = [[0,0,0],[0,0,0],[0,0,0]]
  w = [0.25,0.50,0.25] # Unidimensional weights.
  for i in range(3):
    for j in range(3):
      W[i][j] = w[i]*w[j]
  assert abs(dot(ones, ones, W) - 1) < 1.0e-13;
  return W
  # ----------------------------------------------------------------------

def build_smooth_basis():
  # Returns the basis {B} of the smooth operators, as a list of {n} 3x3
  # arrays, and the corresponding weight window {W}. Also checks that
  # they are orthogonal.
  # 
  # The /basis/ is a list {B[0..n-1]}} of {n} windows 
  # that should be orthogonal under the dotproduct {<|>}.

  W = smooth_weight_window()

  assert n == 9
  B = ( \
    ((00,+1,00),(+1,+1,+1),(00,+1,00)),  # The constant "1" window, i.e. {ones}.
    ((00,00,00),(-1,00,+1),(00,00,00)),  # {X*(1-Y^2)}.
    ((00,-1,00),(00,00,00),(00,+1,00)),  # {Y*(1-X^2)}.
    
    ((00,-1,00),(+2,-2,+2),(00,-1,00)),  # {U = (2*X^2-1)*(1-Y^2)}.
    ((+1,00,-1),(00,00,00),(-1,00,+1)),  # {T = X*Y}.
    ((00,+2,00),(-1,-2,-1),(00,+2,00)),  # {V = (2*Y^2-1)*(1-X^2)}.

    ((+1,00,+1),(00,00,00),(-1,00,-1)),  # {X^2*Y}.
    ((-1,00,+1),(00,00,00),(-1,00,+1)),  # {X*Y^2}.
    ((+1,00,+1),(00,00,00),(+1,00,+1)),  # {X^2*Y^2}.
  )
  
  # Check orthogonality:
  for r in range(n):
    for s in range(n):
      Mrs = dot(B[r], B[s], W)
      if r == s:
        if Mrs <= 0:
          err.write(f"** dot <B{r}|B{r}> = {Mrs}\n"); assert False
      else:
        if Mrs != 0:
          err.write(f"** dot <B{r}|B{s}> = {Mrs}\n"); assert False
  return B, W
  # ----------------------------------------------------------------------

def compute_smooth_ops(Z, B, W):
  # Computes all the 3x3 smoothed operators from the window {Z[0..2][0..2]}.
  # Returns them as a dict.
  
  C = decomp_window(Z, B, W)
  OP = {}

  # Differential operators: computes the estimates of the differential
  # operators at the center of the window, assuming the window samples
  # are the values of a function 
  # 
  #   { J(x,y) = 
  #       C[0] + 
  #       C[1]*x + C[2]*y + 
  #       C[3]*u(x) + C[4]*x*y + C[5]*v(y) + 
  #       C[6]*sx(x,y) + C[7]*sy(x,y) +
  #       C[8]*h(x,y)
  #   }
  #
  # where 
  #
  #   {u(x) = 2*x^2 - 1}
  #   {v(y) = 2*y^2 - 1}
  #   {r(x,y) = x*(x^2 - 2*y^2)}
  #   {s(x,y) = y*(y^2 - 2*x^2)}
  #   {h(x,y) = u(x)*v(y) = (2*x^2 - 1)*(2*y^2 - 1)
  #
  # and the {C} coefficients are computed by {decomp_window}.
  #

  f = C[0]-C[3]-C[5];                   OP["f"] = f
  fx = C[1];                            OP["fx"] = fx
  fy = C[2];                            OP["fy"] = fy
  fxx = 4*C[3];                         OP["fxx"] = fxx
  fxy = C[4];                           OP["fxy"] = fxy
  fyy = 4*C[5];                         OP["fyy"] = fyy
  fxyy = -4*C[6];                        OP["fxyy"] = fxyy
  fxxy = -4*C[7];                        OP["fxxy"] = fxxy
  fxxyy = 16*C[8];                      OP["fxxyy"] = fxxyy

  laplacian = fxx+fyy;                  OP["laplacian"] = laplacian
  orthicity = fxx-fyy;                  OP["orthicity"] = orthicity
  elongation = hypot(2*fxy, orthicity); OP["elongation"] = elongation
  saddleness = 2*fxx*fyy - fxy*fxy;     OP["saddleness"] = saddleness

  # err.write(f"  C[8] = {C[8]}\n")

  # Computes the best weighted least squares fit of a zero-degree
  # polynomial (a constant) to the values in {Z[0..2][0..2]}. Returns
  # that constant as "const1".
  const1 = C[0];  OP["const1"] = const1
  
  # Computes the best weighted least squares fit of 
  # a linear formula {lin1 + linX*X + linY*Y} to the values in {Z[0..2][0..2]}.
  # The coefficient {lin1} is just the weighted average {C[0]}, and the coefficients
  # {linX,linY} are just {fx = C[1],fy = C[2]}. The root mean square of the residual of 
  # that fit is returned as operator "lindev".
  
  lin1 = C[0]; OP["lin1"] = lin1
  linX = C[1]; OP["linX"] = linX
  linY = C[2]; OP["linY"] = linY

  # Computes the best weighted least squares fit of a quadratic formula
  # {quad1 + quadX*X + quadY*Y + quadX2*X^2 + quadXY*X*Y + quadY2*Y^2} to the values in
  # {Z[0..2][0..2]}.
  #
  # Basis elementx {B[0..5]} span the space of all such quadratics, and
  # the other elements {B[5..8]} are orthogonal to them. So the best
  # quadratic fit is the combination of {B[0..5]} with the coefficients
  # {C[0..5]}.  We only need to convert the coefficients {C[0..5]} to
  # the coefficients of the quadratic above.
  #
  # Considering that {B[3] = U = 2*X^2-1} and {B5 = V = 2*Y^2-1},
  # we get that {SUM{ C[r]*B[r] : r \in 0..5 }} is 
  # {(C[0]-C[3]-C[5]) +  C[1]*X + C[2]*Y + 2*C[3]*X^2 + C[4]*X*Y + 2*C[5]*Y^2}
  # Therefore
  
  quad1 = C[0]-C[3]-C[5];  OP["quad1"] = quad1   # Same as the "f" oeprator
  quadX = C[1];            OP["quadX"] = quadX   # Same as "fx".
  quadY = C[2];            OP["quadY"] = quadY   # Same as "fy".
  quadX2 = 2*C[3];         OP["quadX2"] = quadX2 # same as "fxx"/2.
  quadXY = C[4];           OP["quadXY"] = quadXY # same as "fxy".
  quadY2 = 2*C[5];         OP["quadY2"] = quadY2 # same as "fyy"/2.

  return OP
  # ----------------------------------------------------------------------

def compute_avg_dev_ops(Z, W, OPin):
  # Computes the "average" and "deviation" operators.
  # 
  # Also computes the "constdev", "lindev", and "quaddev" deviation
  # operators for the constant, linear, and quadratic approximation
  # coefficients "const1", "lin1", "linX", ... "quadY2" as specified in
  # the dict {OPin}.
  #
  # All these operators are defined to always use the smoothed weight
  # mask, even if the polynomial approximations were computed with the
  # unsmoothed formulas. Thus this procedure assumes that {W} is the
  # basis and weight mask for the smoothed operators.
  
  OP = {}

  # Compute the total weight in the window:
  sum_W = 0
  for i in range(3):
    for j in range(3):
      sum_W += W[i][j]

  # The "average" operator is the weighted average of all pixels
  # in the window:
  ones = [[1,1,1],[1,1,1],[1,1,1]]
  average = dot(ones,Z,W)/sum_W;   OP["average"] = average
  
  # The "deviation", "constdev", "lindev", and "quaddev" are the weighted RMS difference
  # between the nine pixels and the average, constant, linear, and quadratic fits:
  for k in range(4):
    sum_WD2 = 0
    for i in range(3):
      for j in range(3):
        Xij = j - 1
        Yij = i - 1
        if k == 0:
          Vij = average
        elif k == 1:
          Vij = OPin["const1"]
        elif k == 2:
          Vij = \
            OPin["lin1"] + \
            OPin["linX"]*Xij + \
            OPin["linY"]*Yij
        elif k == 3:
          Vij = \
            OPin["quad1"] + \
            OPin["quadX"]*Xij + \
            OPin["quadY"]*Yij + \
            OPin["quadX2"]*Xij*Xij + \
            OPin["quadXY"]*Xij*Yij + \
            OPin["quadY2"]*Yij*Yij
        else:
          assert False
        Dij = Z[i][j] - Vij
        sum_WD2 += W[i][j]*Dij*Dij
    dev = sqrt(sum_WD2/sum_W)
    if k == 0:
      OP["deviation"] = dev
    elif k == 1:
      OP["constdev"] = dev
    elif k == 2:
      OP["lindev"] = dev
    elif k == 3:
      OP["quaddev"] = dev
    else:
      assert False
  return OP
  #----------------------------------------------------------------------

def decomp_window(Z, B, W):
  # Finds linear comb of the the window basis {B[0..n-1]} (assumed orthogonal)
  # that equals the window values {Z[0..2][0..2]}.
  # Returns the coefficients {C[0..n-1]} of the linear combination.
  
  # Compute the values {MD[r] = <B[í]|B[r]>} for {r (assumed diagonal)
  # and indep term {b[r] = <Z|B[r]>}, for {r} in {0..n-1}:
  C = [0]*n
  for r in range(n):
    Mrr = dot(B[r], B[r], W)
    br = dot(Z, B[r], W)
    C[r] = br/Mrr
  return C
  # ----------------------------------------------------------------------

def dot(G, H, W):
  # Given 3x3 windows {G,H}, and a weight window {W}, computes their
  # /dot product/, denoted by {<G|H>}, defined as the sum of
  # {W[i][j]*G[i][j]*H[i][j]} for all sample indices {i,j} in {0..2}.
  # Assumed the weights add to 1.
  
  # If {W} is {None}, assumes all 1s.
  
  sum = 0;
  for i in range(3):
    for j in range(3):
      Wij = 1 if W == None else W[i][j]
      sum += Wij*G[i][j]*H[i][j]
  return sum
  # ----------------------------------------------------------------------

def check_residual(B, W, C, k, val):
  # Checks whether the length of the vector {SUM{C[r]*B[r] : r \in k..n-1}} 
  # is equal to {val}, apart from roundoff.
  sum_L2 = 0
  for r in range(k,n):
    sum_L2 += C[r]*C[r]*dot(B[r], B[r], W)
  if abs(sqrt(sum_L2) - val) > 1.0e-12:
    err.write(f"** residual does not match {sum_L2} {val}\n")
    assert False
  # ----------------------------------------------------------------------

def print_ops(OP):
  for name in sorted(OP.keys()):
    err.write("%-10s %.3f\n" % (name, OP[name]));
  # ----------------------------------------------------------------------

def print_ops4(OPu,OPue,OPs,OPse):
  keys = ( \
    "f", "fx", "fy", "fxx", "fxy", "fyy", "-",
    "fxyy", "fxxy", "fxxyy", "-",
    "laplacian", "orthicity", "elongation", "saddleness", "-",
    "const1", "constdev", "-",
    "lin1", "linX", "linY", "lindev", "-",
    "quad1", "quadX", "quadY", "quadX2", "quadXY", "quadY2", "quaddev", "-",
    "average", "deviation", 
  )
  nk = 0;
  for name in keys:
    if name != "-":
      err.write("%-10s" % name);
      print_ops2(OPu[name], OPue[name])
      print_ops2(OPs[name], OPse[name])
      nk += 1
    err.write("\n")
  assert len(OPu.keys()) == nk
  # ----------------------------------------------------------------------
  
def print_ops2(opv, opve):
  err.write(" %+10.4f" % opv);
  flag = ""
  if opve == None:
    err.write(" %10s" % "??")
  else:
    err.write(" %+10.4f" % opve);
    if abs(opv - opve) > 1.0e-8:
      flag = "<<"
  err.write(" %2s  " % flag)
  # ----------------------------------------------------------------------

def accum(Min, ZMin, Max, ZMax, name, val, Z):
  accum_one(Min, ZMin, -1, name, val, Z);
  accum_one(Max, ZMax, +1, name, val, Z);
  # ----------------------------------------------------------------------

def accum_one(Ext, ZExt, sgn, name, val, Z):
  # The {Ext,ZExt} is asumed to be min if {sgn} is {-1}, max if {sgn} is {+1}.
  if not name in Ext: Ext[name] = -sgn*inf;
  if sgn*val > sgn*Ext[name]:
    Ext[name] = val; ZExt[name] = copy_window(Z)
  elif val == Ext[name]:
    FZ = fluff(Z)
    FZExt = fluff(ZExt[name])
    if FZ < FZExt or (FZ == FZExt and dot(Z,Z,None) < dot(ZExt[name],ZExt[name],None)):
      ZExt[name] = copy_window(Z)
  # ----------------------------------------------------------------------
  
def print_ranges(Min, ZMin, Max, ZMax):
  err.write("=== ranges ===\n")
  for name in Min.keys():
    err.write("%-10s Min = %+8.3f  Max = %+8.3f\n" % (name, Min[name], Max[name]));
    pwin2("ZMin", ZMin[name], "ZMax", ZMax[name])
    err.write("\n");
  # ----------------------------------------------------------------------
  
def fluff(Z): 
  # A measure of how ugly are the numbers in {Z}.
  sum_F2 = 0
  for i in range(3):
    for j in range(3):
      Zij = Z[i][j]
      Fij = Zij - floor(Zij + 1.0e-11)
      sum_F2 += Fij*Fij
  return sqrt(sum_F2/9)
  # ----------------------------------------------------------------------

def pwin(name, Z):
  # Prints the window {Z}. Note that the {Y} axis is down.
  for i in range(3):
    prow(name, Z, i)
    err.write("\n")
  # ----------------------------------------------------------------------
  
def pwin2(name1, Z1, name2, Z2):
  # Prints the windows {Z1,Z2} side by side. Note that the {Y} axis is down.
  for i in range(3):
    prow(name1, Z1, i)
    prow(name2, Z2, i)
    err.write("\n")
  # ----------------------------------------------------------------------
  
def prow(name, Z, i):
  na = name + " = " if i == 1 else (" "* len(name)) + "   "
  err.write("  " + na + "[")
  for j in range(3):
    s = ("%+6.3f" % Z[i][j])
    s = re.sub(r'[+]0[.]', " 0.", s)
    k = re.search(r'[.]?[0]*$', s).start()
    n = len(s)
    s = s[0:k]
    if j < 2: s += ","
    s = (s + " "*(n-k))
    err.write(s)
  err.write("]")
  # ----------------------------------------------------------------------

def copy_window(Z):
  ZC = [[0,0,0],[0,0,0],[0,0,0]]
  for i in range(3):
    for j in range(3):
      ZC[i][j] = Z[i][j]
  return ZC
  # ----------------------------------------------------------------------

test_basis(False)
test_basis(True)
test_ops()
test_enum_windows(False)
test_enum_windows(True)
