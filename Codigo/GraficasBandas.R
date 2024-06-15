#Some libraries that are required to plot
library(plot3D)
library(pals)

#We set up the graph as to see from every angle (this is just for me to check)
par(mfrow=c(1,2), mar = c(1.5, 3.2, 4.5, 1) + 0.1, xpd=T)

#We need to color the surfaces using the same scale, so I'll write a function 
#for a reasonable color assignment
colors = function(k) {
  n = eigenplot*30
  paleta = jet(n)
  max = max(E[,,eigenplot])
  min = min(E[,,1])
  localmax = max(E[,,k])
  localmin = min(E[,,k])
  cota1 = floor(abs(n*localmin/max))
  cota2 = floor(abs(n*localmax/max))
  return(paleta[min(c(cota1,cota2)):max(c(cota1,cota2))])
}


#Some parameters that control the look of the graphs
shadows = 0.2 #[0,1]
transparency = 0.7 #[0,1]

#Now we plot the energy surfaces
persp3D(z = E[,,eigenplot], x = K$x, y = K$y, col = colors(eigenplot), shade = shadows,
        alpha=transparency, facets = T, scale = T, zlim=c(-1,max(E)),
        phi=30, theta=120, colkey = list(plot = FALSE),border="black",lwd=0.1,
        xlab="",ylab="",zlab="", ticktype="detailed",cex.axis=0.8)

for(k in (eigenplot-1):1) persp3D(z = E[,,k], x = K$x, y = K$y,
                                  col = colors(k), shade = shadows, alpha=transparency,
                                  border="black",lwd=0.1, facets = T, scale = T,
                                  add=T, colkey = list(plot = FALSE))


#We add some necessary labels
text3D(0.7, 1.3, 0, expression(K[x]*a/pi), add=T, cex=0.85)
text3D(1.65, -0.1, 0, expression(K[y]*a/pi), add=T, cex=0.85)
text3D(2.1, 0.3, 5.75, expression(E (eV)), add=T, cex=0.85)
title("Superficies de energía", cex.main=1.5, outer=T, line=-6.5)
title("Celdilla unidad (Kronig-Penney 2D)", cex.main=0.9,
      outer=T, line=-7.5, font.main=1)

#Now we plot the energy surfaces
persp3D(z = E[,,eigenplot], x = K$x, y = K$y, col = colors(eigenplot), shade = shadows,
        alpha=transparency, facets = T, scale = T, zlim=c(-1,max(E)), border="black", lwd=0.05,
        phi=0, theta=0, colkey = list(plot = FALSE),
        xlab="",ylab="",zlab="",cex.axis=0.8)

#And again some labels
text3D(-0.1,-2.3,0, expression(K[x]*a/pi), add=T, cex=0.85)
text3D(-0.76,-2.2,0.4, expression(E (eV)), 
       outer=T, add=T, cex=0.8)

#The rest of the plot
for(k in (eigenplot-1):1) persp3D(z = E[,,k], x = K$x, y = K$y,
                                  border="black",lwd=0.05,
                                  col = colors(k), shade = shadows, alpha=transparency,
                                  facets = T, scale = T, add=T, colkey = list(plot = FALSE))

#The color key goes in the middle, for example
colkey(col = jet(100),
       clim = c(min(E[,,1]), max(E[,,eigenplot])), 
       at = c(0:3), add = T, side=2, cex.axis=0.8, length=0.3)