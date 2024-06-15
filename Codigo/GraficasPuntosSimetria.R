#Now we will plot the trajectory gamma -> X -> M -> gamma.
#Load the LaTeX package
library(latex2exp)

#Make an adequate grid for the graphs
par(mfrow=c(1,2), xpd=T, mar=c(11,3,7,1), bg="white")

#Generate some colors
cols = c("#DB4437","#4285F4", "#0F9D58")

#We plot the energy surfaces from above and show the path
image(E[,,eigenplot],col=jet(101),asp=1,xaxt="n",yaxt="n")
arrows(1,1,0.51,0.51,add=T,lwd=2,length=0.08,col="white") #Path 1, X to M
arrows(0.5,0.5,1,0.5,add=T,lwd=2,length=0.08,col="white") #Path 2, X to M
arrows(1,0.51,1,1,add=T,lwd=2,length=0.08,col="white") #Path 3, X to M
text(0.6,0.7,labels=TeX("$M\\rightarrow \\Gamma$"),
     col="black",add=T,cex=1)
text(0.75,0.45,labels=TeX("$\\Gamma\\rightarrow X$"),
     col="black",add=T,cex=1)
text(0.9,0.65,labels=TeX("$X\\rightarrow M$"),
     col="black",add=T,cex=1)



#Extract the size of the matrix
size = dim(E[,,1])[1]/2

#Create a blank plot
plot(0,0,type="n", ylim=c(min(E),max(E)), xlim=c(0,size*3),
     xlab="", xaxt="n")
mtext("E (eV)", side=2, line=2)

#And plot the profiles
for(k in eigenplot:1) {
  #We create a subset of the relevant energy surface
  E_graph = E[size:(size*2),size:(size*2),k]
  #We extract each path
  path1 = E_graph[1,]
  path2 = E_graph[,size]
  path3 = rev(diag(E_graph))
  #Append them all
  path_total = c(path1,path2,path3)
  #And add them to the plot
  lines(path_total,lwd=3,lty="solid",col=cols[k])
}

#We add the lines showing where each path ends
text(size/2,0.5,TeX("$\\Gamma\\rightarrow X$"),bg="white")
text(3/2*size,0.5,TeX("$X \\rightarrow M$"),bg="white")
text(5/2*size,0.5,TeX("$M \\rightarrow \\Gamma$"),bg="white")
abline(v=size+3,lty="dotted",col="gray",lwd=1,xpd=F)
abline(v=2*size+3,lty="dotted",col="gray",lwd=1,xpd=F)
abline(v=3*size+2,lty="dotted",col="gray",lwd=1,xpd=F)

