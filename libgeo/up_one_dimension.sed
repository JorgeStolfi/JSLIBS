#! /bin/sed -f
# Last edited on 2024-08-30 11:39:18 by stolfi


s:3x3:4x4:g
s:r3:r4:g
s:\[3\]:[4]:g
s:R\^3:R\^4:g
s:\]\^3:]^4:g
s:S\^3:S\^4:g
s:3-ball:4-ball:g
s:\bplane\b:hyperplane:g 
s:\bplanes\b:hyperplanes:g 
s:triangle:tetrahedron:g
s:collinear:coplanar:g
s:quadrilateral:cuboid:g

s:2x2:3x3:g
s:r2:r3:g
s:\[2\]:[3]:g
s:R\^2:R\^3:g
s:\]\^2:]^3:g
s:S\^2:S\^3:g
s:2-ball:3-ball:g
s:\bline\b:plane:g
s:\blines\b:planes:g

s:\[1\]:[2]:g
s:R\^1:R\^2:g
s:S\^1:S\^2:g
