
Exercise 1:

Task 1:

Cos(θ) = (x · y) / (||x|| * ||y||)
(x · y) = sqrt(2) * 1 + 1 * 1 + 0 * 1 = sqrt(2) + 1

|| x || = sqrt(sqrt(2)**2 + 1**2 + 0**2) = sqrt(3)
|| y || = sqrt(1**2 + 1**2 + 1**2) = sqrt(3)
(||x|| * ||y||) = sqrt(3) * sqrt(3) = 3
Cos(θ) = (sqrt(2) + 1 )/ 3 = 0.80473785412

Task 2:
Use the cross product to find a vector perpendicular to x and y
z = x × y 
[
    x2 * y3 - x3 * y2
    x3 * y1 - x1 * y3
    x1 * y2 - x2 * y1
]
= 
[    
    1 * 1 - 0 * 1
    0 * 1 - sqrt(2) * 1
    sqrt(2) * 1 - 1 * 1
]
=
[1,-sqrt(2), sqrt(2) - 1] = [1, -1.41421356, 0.41421356]

To normalise z we divid each element of z by its magnitude
|| z || = sqrt(1**2 + (-sqrt(2))**2 + (sqrt(2) - 1) ** 2)

 = sqrt(6- 2 * sqrt(2)) = 1.78089103408

Now we compute z / || z ||
[1,-sqrt(2), sqrt(2) - 1] / sqrt(6- 2 * sqrt(2)) = [ 0.56151667 -0.79410449  0.23258782]


Task 3:
A = [
    1 1 1
    2 2 1
    -1 -3 -3
]
z = [ 0.56151667 -0.79410449  0.23258782]

u = Az = [0, -0.2326, 1.1230]






