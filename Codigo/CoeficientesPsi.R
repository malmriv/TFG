#El índice de vectors[,i] indica el índice de la función de onda.
#Comprobamos que a energías bajas dominan los coeficientes pequeños
library(viridis)
library(latex2exp)
library(fields)
par(mar=c(4.2,4.2,4.2,4.2))

image2 = function(mat1) {
  n = length(mat1[1,])
  mat1 <- apply(mat1, 2, rev)
  image.plot(1:n,1:n,t(mat1),col=viridis(100),yaxt="n",xaxt="n",asp=1)
  mtext(side=1,TeX("Índice de la banda"),line=2.1,cex=1.2)
  mtext(side=2,TeX("$|b(n_x,n_y)|^2$"),line=1.5,cex=1.2)
  axis(side=1, at=seq(0,sqrt(length(vectors)),by=10), labels = T)
}
coeffmatrix = abs(vectors)^2
image2(coeffmatrix)


options(digits=25)
for(j in 1:length(coeffmatrix[1,])) print(sum(coeffmatrix[,j]))

par(mfrow=c(2,3))
for(k in 1:eigenplot) image(as.matrix(blur(as.im(wavefunction[,,k]),sigma=0.5)),asp=1,col=jet(150),xlab="x/a",ylab="y/a")
  