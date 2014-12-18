### Script for Phylogenetic Analysis of Pitches in Anura ###
### Authors: David Lucas Röhr; Gustavo B. Paterno; Felipe Camurugi; Flora A. Juncá ; Adrian A. Garda
### data created: 25/10/13
### Updated: 09/01/14
### by: Gustavo Paterno
library(ape);library(caper);library(phytools)
library(nlme);library(dplyr);library(ggplot2);library(knitr)
library(picante)
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


### DATA ###
mat <- read.table("raw.data_Rhor_JEB.txt",h=T)
str(mat)

### Análise exploratória:
summarise(group_by(mat,environment),meanDF=mean(DF),seDF=stderr(DF),
                meanSVL=mean(SVL),seSVL=stderr(SVL))
summary(mat)
summary(subset(mat,mat$amb=="len"))
summary(subset(mat,mat$amb=="lot"))
sd(mat$crc)
sd(mat$fd)
med.fd <- with(mat,tapply(fd,amb,mean))
err.fd <- with(mat,tapply(fd,amb,sd))
med.crc <- with(mat,tapply(crc,amb,mean))
err.crc <- with(mat,tapply(crc,amb,sd))

med.fd[2] - med.fd[1]  ### diferença da FD entre ambientes

# Subset só com dados Len / Lot
mat.len <- subset(mat,mat$environment=="still")
mat.lot <- subset(mat,mat$environment=="running")
mat$sp
### Shapiro test for lentic habitat:
shapiro.test(mat.len$DF)
shapiro.test(mat.len$logDF)
shapiro.test(mat.len$SVL)
shapiro.test(mat.len$logSVL)

### Shapiro test for lotic habitat:
shapiro.test(mat.lot$DF)
shapiro.test(mat.lot$logDF)
shapiro.test(mat.lot$SVL)
shapiro.test(mat.lot$logSVL)

### Data: Phylogenetic Tree (Pyron, 2011)
arvore<-read.tree("amph_2014.tre")

anura<-drop.tip(arvore,as.character(mat[,1]))               ### droping species from tree
filo<-(drop.tip(arvore,anura$tip.label))
filo$node.label ### Verificando node labels
filo$node.label <- makeLabel(filo)$node.label ### gerando node.labels únicos
plot(filo,cex=.1)
length(filo$tip.label)
write.tree(filo,"phy_Rohr_JEB_2014.tre")
write.tree(arvore,"amph_2014.tre")

anura_tree<-read.tree("anura_tree.tre")
write.table(filo$tip.label,"species_list.txt",row.names=F)

### Checking for absent species in data:
sort(filo$tip.label) == sort(mat$sp)      ### Quais são os erros de digitação
sum(sort(filo$tip.label) != sort(mat$sp)) ### Quantos são os erros de digitação

### Data preparation for pgls:
comp.data <- comparative.data(phy=filo,data=mat,names.col="sp",vcv=T,vcv.dim=3)
comp.data.runn <- comparative.data(phy=filo,data=subset(mat,environment=="running"),names.col="sp",vcv=T,vcv.dim=3)
comp.data.still <- comparative.data(phy=filo,data=subset(mat,environment=="still"),names.col="sp",vcv=T,vcv.dim=3)

######################################### Analysis ########################################################
### Phylogenetic Signal:
comp.data$data$phyres <- mod.pgls$phyres
k.signal <- multiPhylosignal(select(comp.data$data,logDF,logSVL,phyres),comp.data$phy,reps=999)
kable(k.signal)


### PGLS with with lambda adjusted by maximum likelihood 
mod.pgls <- pgls(logDF ~ environment*logSVL, data=comp.data,lambda=1)

### Anova Table for each model:
mod.pgls1.ano <- anova(mod.pgls1)
mod.pgls2.ano <- anova(mod.pgls2)
mod.pgls3.ano <- anova(mod.pgls3)
anova(mod.pgls1)
summary(mod.pgls1)
summary(mod.pgls2)
summary(mod.pgls3)

pgls1.cof <- coef(mod.pgls1)
pgls2.cof <- coef(mod.pgls2)
pgls3.cof <- coef(mod.pgls3)

# Best model:
mod.pgls2.ano

### Coeficientes das retas de regressão
len.cof <- pgls3.cof[c(1,2)]
lot.cof <- c(c(pgls3.cof[1]+pgls3.cof[3]),pgls3.cof[2])
# Diferença entre os interceptos (len vs lot) em unidade original (Hertz)
diff.intercep <- exp(lot.cof[1]) - exp(len.cof[1])

################################ Analysis with Family ########################################################
table(mat$fam) > 30
big.fam <- c("Bufonidae","Centrolenidae","Hylidae","Hyperoliidae",
             "Leptodactylidae","Microhylidae","Ranidae")
mat.fam <- droplevels(filter(mat,
                  fam == big.fam[1]|
                  fam == big.fam[3]))
                  

str(mat.fam)
table(mat.fam$fam)
### Preparaing data for pgls:
c.dat.fam <- comparative.data(phy=filo,data=mat.fam,names.col="sp",vcv=T,vcv.dim=3)
mod.fam <- pgls(logfd ~ amb*logcrc+fam*amb, data=c.dat.fam,lambda="ML")
anova(mod.fam)
summary(mod.fam)
g.fam <- ggplot(mat.fam,aes(x=logcrc,y=logfd,colour=fam:amb,group=amb))+
          facet_grid(~fam)+
          geom_point(alpha=.6)+
          geom_smooth(method="lm",se=F);g.fam

g.fam2 <- ggplot(mat,aes(x=logcrc,y=logfd,colour=fam:amb,group=fam))+
          #facet_grid(~fam)+
          geom_point(alpha=.6)+
          geom_smooth(method="lm",se=F);g.fam2

### PGLS family:
mat.buf <- droplevels(filter(mat.fam,fam==big.fam[1]))
mat.hyl <- droplevels(filter(mat.fam,fam==big.fam[3]))
mat.ran <- droplevels(filter(mat.fam,fam==big.fam[7]))

c.dat.buf <- comparative.data(phy=filo,data=mat.buf,names.col="sp",vcv=T,vcv.dim=3)
c.dat.hyl <- comparative.data(phy=filo,data=mat.hyl,names.col="sp",vcv=T,vcv.dim=3)
c.dat.ran <- comparative.data(phy=filo,data=mat.ran,names.col="sp",vcv=T,vcv.dim=3)

mod.buf <- pgls(logfd ~ amb*logcrc, data=c.dat.buf,lambda="ML")
mod.hyl <- pgls(logfd ~ amb*logcrc, data=c.dat.hyl,lambda="ML")
mod.ran <- pgls(logfd ~ amb*logcrc, data=c.dat.ran,lambda="ML")

anova(mod.buf)
anova(mod.hyl)
anova(mod.ran)

### Martins, E. P. and Hansen, T. F. (1997) 
mat <- arrange(mat,sp)
rownames(mat) <- mat$sp
mydata <- mat[match(filo$tip.label,rownames(mat)),]

ou.anura<-corMartins(1,phy=filo,fixed=T) 
ou.gls<-gls(logfd ~ amb*logcrc,correlation=ou.anura,data=mydata) 
summary(ou.gls)
anova(ou.gls)
sum(mat$sp == anura_data1$phy$tip.label)


