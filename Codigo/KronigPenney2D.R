# Title: 2D Kronig-Penney numerical simulation, as explained by 
# Pavelich and Marsiglio (2016). Preliminary work for my TFG.
# Author: Manuel Almagro Rivas (malmriv@correo.ugr.es)

#Load necessary libraries
library(plot3D)

# 2D finite well with open boundary conditions (Bloch
# conditions will be applied later)
periodicwell2D = function(N, v, p1, p2) {
  
  # A more convenient way of navigating the quantum number list:
  # converts {1,2,3,4,...} -> {0,1,-1,2,-2,...}
  j = 1:(N + 1)
  a = (1 + (2 * j - 1) * (-1) ^ j) / 4
  
  # We generate pairwise combinations
  n_matrix = expand.grid(a, a)
  n = as.matrix(n_matrix)
  nx = n[, 1]
  ny = n[, 2]
  
  # We are interested in ordering the quantum number list so that
  # we begin by computing the lower energy eigenstates
  n2 = rep(0, (N + 1) * (N + 1))
  n = cbind(n2, nx, ny)
  
  # Appending n^2 column
  n[, 1] = n[, 2]^2 + n[, 3]^2
  
  # Sorting by the n^2 column
  n = n[order(n[, 1]), ]
  
  # Sorting by n^2 to give increasing energy encoding
  I = order(n[, 1])
  T_seq = n[I, 1]
  n = n[I, ]
  nx = n[, 2]
  ny = n[, 3]
  n2 = n[, 1]
  
  # Populating diagonal matrix elements
  p = p2 - p1
  h = diag(4 * n2 + v * p^2)
  
  # Precomputing off-diagonal integrals
  off = rep(0, N + 1)
  for (m in 1:(2 * N + 1)) {
    off[m] = 1 / (2 * pi * 1i * (m - N - 1)) * (-exp(1i * 2 * pi * (m - N - 1) * p1) + exp(1i * 2 * pi * (m - N - 1) * p2))
  }
  
  ipart = 0
  jpart = 0
  
  # Populating off-diagonal elements
  for (i in 1:((N + 1) * (N + 1) - 1)) {
    for (j in (i + 1):((N + 1) * (N + 1))) {
      # NA's should be avoided!
      if (!is.na(nx[i]) && !is.na(nx[j]) && nx[i] == nx[j]) {
        ipart = p
      } else {
        ipart = off[nx[i] - nx[j] + N + 1]
      }
      
      if (!is.na(ny[i]) && !is.na(ny[j]) && ny[i] == ny[j]) {
        jpart = p
      } else {
        jpart = off[ny[i] - ny[j] + N + 1]
      }
      
      h[i, j] = h[i, j] + v * ipart * jpart
    }
  }
  
  # Fill the diagonals and take the conjugate transpose for the lower diag
  h = h + t(h)
  diag(h) = diag(h) / 2
  
  #Return the dimensionless hamiltonian
  return(h)
}

#Bloch condition is applied to the 2D periodic well
#   The whole unit cell is covered, unlike the code given by
#   Pavelich and Marsiglio, which only covers high-symmetry points.
#
#   N  = number of basis states in 1D (total matrix size is N+1 x N+1)
#   v  = depth of well, should be negative
#   p1 = start of the well in units of a
#   p2 = end of the well in units of a
#   iter = number of points in K-space to sample from
#   E  = eigenvalues of resultant matrix diagonalizations
#   Kx, Ky = array of points along K_x, K_y axis
#   nx, ny = quantum numbers of plane-wave basis states
#   h = dimensionless Hamiltonian

# The parameters are set here
N = 6 #simulation complexity grows as O(N^2)
a = 1.611575e-09
p1 = 0.3088397
p2 = 1-p1
iter = 41 #number of K points to sample in each dimension
#this HAS TO be odd, I didn't want to make the code
#too complicated when graphing and just assumed this
#to be odd
eigenplot = 3 #number of energy surfaces to plot
E_ISW = 4*pi^2*hbar^2/(2*m0*a^2)
conversionfactor = E_ISW/1.609e-19
v = -5.018602/conversionfactor

#Parameters to study (dot to dot distance, abs and rel)
dot2dotrel = 2*p1*100
dot2dotabs = 2*p1*a

# Initial calculation of matrix elements
h = periodicwell2D(N, v, p1, p2)

# This block is just like what we did before, check previous
# comments for more details
j = 1:(N + 1)
b = (1 + (2 * j - 1) * (-1) ^ j) / 4
n_matrix = expand.grid(b, b)
n = as.matrix(n_matrix)
nx = n[, 1]
ny = n[, 2]
n2 = rep(0, (N + 1) * (N + 1))
n = cbind(n2, nx, ny)
n[, 1] = n[, 2]^2 + n[, 3]^2
n = n[order(n[, 1]), ]
I = order(n[, 1])
T_seq = n[I, 1]
n = n[I, ]
nx = n[, 2]
ny = n[, 3]
n2 = n[, 1]

# This is new: we generate a rectangular K-space grid
Kx = seq(-pi, pi, length.out = iter + 1)
Ky = seq(-pi, pi, length.out = iter + 1)
K = mesh(Kx,Ky)

# Matrix indexing parameter, so that every quantum number
# is considered in the calculation
m = 1:((N + 1) * (N + 1))

#For every K value we will have several energy surfaces. We'll compute
#eigenvalues for each (Kx,Ky) point, and then plot them
E = array(0, dim = c(iter+1, iter+1, eigenplot))

# We initialise the progress bar
progressbar = txtProgressBar(min = 0, max = (iter+1), style = 3, width = 50, char = "|")

for(i in 1:(iter+1)) {
  setTxtProgressBar(progressbar, i) #Show progress to user (me)
  for(j in 1:(iter+1)) {
    #Bloch correction to hamiltonian:
    bloch = (4/pi) * (nx[m] * Kx[i] + ny[m] * Ky[j]) +
      (Kx[i]^2 + Ky[j]^2) / (pi^2)
    cEb = eigen(h + diag(bloch))
    Ebloch = sort(cEb$values)
    for(k in eigenplot:1) E[i,j,k] = Ebloch[k]
  }
}

#The energies need to be converted to eV, instead of our escalated energy unit
for(k in 1:eigenplot) E[,,k] = E[,,k]*conversionfactor

#Now that we've used the K-grid, we'll make it non-dimensional too
K$x = K$x/pi
K$y = K$y/pi
