#Some libraries that are required to plot
library(plot3D)
library(pals)
library(spatstat)

#The following variable controls if the plots should be generated (TRUE) or not (FALSE)
shouldPlot = F

#We will plot a unit cell with the well we defined previously.
#Generate a set of points to plot
x = seq(0, a, length.out = iter + 1)
y = seq(0, a, length.out = iter + 1)
unitcell = mesh(x,y)

#Create an array to store the points that make up the well
well = array(0, dim=c(iter+1,iter+1))

#We define the 2D well
V = function(x,y) {
  if(x>(p1*a) && x<(p2*a) && y>(p1*a) && y<(p2*a)) return(v*conversionfactor)
  else{return(0)}
}

#We compute each point
for(i in 1:(iter+1)) {
  for(j in 1:(iter+1)) {
    well[i,j] = V(x[i], y[j])
  }
}

#We'll need the hamiltonian matrix to be diagonalised, arrange the
#eigenvalues by increasing order, and then order the eigenvectors
#following the exact same order. This function does that:
OrderedEigenvectors = function(A) {
  eigenA = eigen(A)$vectors
  B = matrix(0, nrow=nrow(A),ncol=ncol(A))
  order = sort(Re(eigen(A)$values), index.return=T, decreasing=F)$ix
  for(k in 1:ncol(B)) {
    for(l in 1:nrow(B)) B[l,k] = eigenA[l,order[k]]
  } 
  return(B)
}

#Each point in the grid will have as many wavefunctions as
#energy surfaces we want to plot. For one surface, we'll need an array
#with dimensions (iter+1) x (iter+1) x (N+1)^2
coeffs = array(0,dim=c(iter+1,iter+1,(N+1)^2,eigenplot))

for(i in 1:(iter+1)) {
  setTxtProgressBar(progressbar, i)
  for(j in 1:(iter+1)) {
    #Bloch correction to hamiltonian:
    bloch = (4/pi) * (nx[m] * Kx[i] + ny[m] * Ky[j]) +
      (Kx[i]^2 + Ky[j]^2) / (pi^2)
    vectors = OrderedEigenvectors(h + diag(bloch))
    for(k in 1:eigenplot)
    coeffs[i,j,,k] = vectors[,k]
  }
}

#The values stored in psi are actually the coordinates in the 
#plane wave basis we chose (i.e. the coefficients that go with each element)
wavefunction = array(0, dim=c(iter+1,iter+1,eigenplot))

#We compute each point
for(b in 1:eigenplot) {
for(i in 1:(iter+1)) {
  for(j in 1:(iter+1)) {
    xi = x[i]
    yi = y[j]
    for(k in 1:(N+1)^2) {
      wavefunction[i,j,b] = wavefunction[i,j,b] + 
        coeffs[i,j,k,b]*exp(2*pi*1i*n[k,2]*xi/a)*exp(2*pi*1i*n[k,3]*yi/a)/a
    }
  }
}
}

#Compute the square
wavefunction = abs(wavefunction)^2

#PLOTS GO HERE.
if(shouldPlot) {
#Some parameters that control the look of the graphs
par(mfrow=c(1,2), mar = c(1.5, 2.5, 4.5, 1) + 0.1, xpd=T)
shadows = 0 #[0,1]
transparency = 0.05 #[0,1]

#The well plot
persp3D(z = well, x = unitcell$x/a, y = unitcell$y/a, col = rgb(0,0.7,0.6), shade = shadows,
        alpha=transparency, facets = T, scale = T, zlim=c(v,0), border="black",
        lwd=0.005, phi=30, theta=50, colkey = list(plot = FALSE),
        xlab="",ylab="",zlab="",cex.axis=0.75, ticktype="detailed")

#We add some necessary labels
text3D(1.2, 1/2, v, expression(x/a), add=T, cex=0.85)
text3D(0.5, -0.3, v, expression(y/a), add=T, cex=0.85)
text3D(-0.25, -0.1, 0.2, expression(V (eV)), add=T, cex=0.85)


#The wavefunction plot
persp3D(z = wavefunction[,,1], x = unitcell$x/a, y = unitcell$y/a, col = jet(100), shade = 0,
        alpha=1, facets = T, scale = T, zlim=c(0,max(wavefunction[,,1])),
        phi=30, theta=50, colkey = list(plot = FALSE),border="black",lwd=0.05,
        xlab="",ylab="",zlab="", cex.axis=0.75, zaxt="n")

#We add some necessary labels
text3D(1.2, 1/2, 0, expression(x/a), add=T, cex=0.85)
text3D(0.5, -0.3, 0, expression(y/a), add=T, cex=0.85)
text3D(0.2, -0.2, max(wavefunction), expression(abs(psi)^2), add=T, cex=1)

#Another visualization for the wavefunction
#dev.new()
#wavefunction_blur = as.matrix(blur(as.im(wavefunction[,,1]),sigma=0.5))
#image(wavefunction_blur,asp=1,col=jet(150))
}