# transcription level across all the genes
# Xiaoji Sun
# 12/11/13

# Read wiggle files ######################################################################################
t1 = read.table('mRNA_treplication.coverage.F.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)
t2 = read.table('mRNA_treplication.coverage.Rpos.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)

t1 = read.table('mRNA_tsc_recomb.coverage.F.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)
t2 = read.table('mRNA_tsc_recomb.coverage.Rpos.rpkm.wig', skip=1, fill=T, stringsAsFactors=F)

index1 = which(t1[,3]!="")
index2 = which(t2[,3]!="")

data1 = list()
for (i in 1:16)
{
  data1[[i]] = t1[(index1[i]+1):(index1[i+1]-1),1:2]
}

data2 = list()
for (i in 1:16)
{
  data2[[i]] = t2[(index2[i]+1):(index2[i+1]-1),1:2]
}

# load the gff file (saved as txt by excel)
out=read.table('~/Documents/LAB/data_analysis/sk1_genome/sk1_sgrp_genes.gff')
# extract the GENEs + 650bp on either side
newmatrix=matrix(NA,ncol=5,nrow=1)
for (i in 1:(nrow(out)-1))
{
  if (strsplit(as.character(out[i,5]),'')[[1]][1]=='Y')
  {
    print (i)
    newline=out[i,]
    newmatrix=rbind(newmatrix,newline)
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
mid=newmatrix

midmatrix=as.data.frame(matrix(NA,ncol=4,nrow=1))
for (i in 1:nrow(mid))
{
  newline=c(mid[i,1],as.numeric(mid[i,2])-650,as.numeric(mid[i,3])+650,mid[i,4])
  midmatrix=rbind(midmatrix,newline)
}
midmatrix=midmatrix[2:nrow(midmatrix),]

# ORF regions
midmatrix=cbind(midmatrix, 2000/(as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])))

# STRAND==+
midmatrix_w=midmatrix[which(midmatrix[,4]=='+'),]

# STRAND==-
midmatrix_c=midmatrix[which(midmatrix[,4]=='-'),]

# normalize genes on the watson strand
red1=data1
for (i in 1:16)
{
  red1[[i]][,1]=as.numeric(red1[[i]][,1])
  red1[[i]][,2]=as.numeric(red1[[i]][,2])
}
alldata_mid_w=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  ## strand
  chr=midmatrix_w[which(midmatrix_w[,1]==c),]
  chr[,2]=as.numeric(chr[,2])
  chr[,3]=as.numeric(chr[,3])
  if (nrow(chr)!=0)
  {wigmatrix=matrix(NA,ncol=2,nrow=1)
   for (i in 1:nrow(chr))
   {
     start=as.numeric(chr[i,2])
     end=as.numeric(chr[i,3])
     if ((length(which(as.numeric(red1[[k]][,1])>=start & as.numeric(red1[[k]][,1])<=end)))!=0)
     {
       tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,5]
       tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
       tmp=cbind(tmp1,tmp2)
       wigmatrix=rbind(wigmatrix,tmp)
     }
   }
   wigmatrix=wigmatrix[2:nrow(wigmatrix),]
   
   # mean
   wig_sort=wigmatrix[order(wigmatrix[,1]),]
   wig_sort[,1]=round(wig_sort[,1])
   mean_matrix_w=as.data.frame(matrix(NA,ncol=2,nrow=2001))
   for (i in 0:2000)
   {
     mean_matrix_w[i+1,1]=i
     mean_matrix_w[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
   }
   #plot(mean_matrix, pch=16,cex=0.6,col='blue')
   
   # smooth
   bp=5
   sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
   for (i in seq(1,2000, bp))
   {
     sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
     sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
   }
   alldata_mid_w[[k]]=sm_matrix
   #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
  }
}

# normalize genes on the crick strand
alldata_mid_c=list()
red1=data2
for (i in 1:16)
{
  red1[[i]][,1]=as.numeric(red1[[i]][,1])
  red1[[i]][,2]=as.numeric(red1[[i]][,2])
}
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  ## strand
  chr=midmatrix_c[which(midmatrix_c[,1]==c),]
  if (nrow(chr)!=0)
  {wigmatrix=matrix(NA,ncol=2,nrow=1)
   for (i in 1:nrow(chr))
   {
     start=as.numeric(chr[i,2])
     end=as.numeric(chr[i,3])
     if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
     {
       tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,5]
       tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
       tmp=cbind(tmp1,tmp2)
       wigmatrix=rbind(wigmatrix,tmp)
     }
   }
   wigmatrix=wigmatrix[2:nrow(wigmatrix),]
   
   # mean
   wig_sort=wigmatrix[order(wigmatrix[,1]),]
   wig_sort[,1]=round(wig_sort[,1])
   mean_matrix_c=as.data.frame(matrix(NA,ncol=2,nrow=2001))
   for (i in 0:2000)
   {
     mean_matrix_c[i+1,1]=i
     mean_matrix_c[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
   }
   #plot(mean_matrix, pch=16,cex=0.6,col='blue')
   
   # smooth
   bp=5
   sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
   for (i in seq(1,2000, bp))
   {
     sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
     sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
   }
   alldata_mid_c[[k]]=sm_matrix
   #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
  }
}



# Combine both strands
alldata_mid=list()
for (k in 1:16)
{
  start=alldata_mid_w[[k]][,1]
  w=alldata_mid_w[[k]][,2]
  c=alldata_mid_c[[k]][nrow(alldata_mid_c[[k]]):1,2]
  temp=as.data.frame(matrix(NA,ncol=2,nrow=1))
  for (i in 1:400)
  {
    newline=c(start[i],(as.numeric(w[i])+as.numeric(c[i]))/2)
    temp=rbind(temp,newline)
  }
  alldata_mid[[k]]=temp[2:nrow(temp),]
}

#
mean_matrix_all=cbind(mean_matrix_w, mean_matrix_c[nrow(mean_matrix_c):1,2])
mean_matrix_all=cbind(mean_matrix_all, rowMeans(mean_matrix_all[,2:3]))
mean_matrix=mean_matrix_all[,c(1,4)]

# plot each chromosome
color_s=c('red','yellow','orange','blue','purple','springgreen','brown','burlywood','pink','grey','skyblue','turquoise','dark green','violetred','magenta','gold')
color=color_s
png("test.png", width=700, height = 500)
plot(alldata_mid[[1]][,1],alldata_mid[[1]][,2], pch=16,col=color[1],
     xlab='Normalized 5\'-3\' ORF (bp)', ylab=Name)
for (i in 2:16)
{
  points(alldata_mid[[i]][,1],alldata_mid[[i]][,2],col=color[i], pch=16)
}
legend(x=0, y=4, col=color, pch=16, legend=1:16, horiz=T, cex=0.85)
dev.off()


# plot average signals
png("test.png", width=1200, height = 1000)
plot(mean_matrix, pch=16, col='blue', xlab='Normalized 5\'-3\' ORF (bp)', ylab=Name, type='l',lwd=2)
dev.off()




###############################################################################################
# if use transcripts ends!
# load the gff file (saved as txt by excel)
out=read.table('~/Documents/LAB/data_analysis/Transcription/TransEnds_genes.txt',header=T)
w2=out[which(out[,4]=='+'),]
c2=out[which(out[,4]=='-'),]
n = 200
w_matrix=data.frame(matrix(0,nrow=nrow(w2),ncol=2*n))
for (i in 1:nrow(w2))
{
  print (i)
  r = as.numeric(w2[i,'chr'])
  #chr_data = red1_data[[r]]
  chr_data = data1[[r]]
  mid = as.numeric(w2[i,'TransEnd'])
  start = as.numeric(w2[i,'TransEnd'])-n
  end = as.numeric(w2[i,'TransEnd'])+n
  tmp = chr_data[which(as.numeric(chr_data[,1])>=start & as.numeric(chr_data[,1])<end),]
  tmp_vector=rep(0,length=2*n)
  for (j in 1:nrow(tmp))
  {
    pos=as.numeric(tmp[j,1])-mid+n+1
    tmp_vector[pos]=as.numeric(tmp[j,2])
  }
  w_matrix[i,]=tmp_vector
}
write.table(w_matrix,'rpkm_200_w.txt', sep='\t', quote=F, col.names=F, row.names=F)

n = 200
c_matrix=data.frame(matrix(0,nrow=nrow(c2),ncol=2*n))
for (i in 1:nrow(c2))
{
  print (i)
  r = as.numeric(c2[i,'chr'])
  #chr_data = red1_data[[r]]
  chr_data = data2[[r]]
  mid = as.numeric(c2[i,'TransEnd'])
  start = as.numeric(c2[i,'TransEnd'])-n
  end = as.numeric(c2[i,'TransEnd'])+n
  tmp = chr_data[which(as.numeric(chr_data[,1])>=start & as.numeric(chr_data[,1])<end),]
  tmp_vector=rep(0,length=2*n)
  for (j in 1:nrow(tmp))
  {
    pos=as.numeric(tmp[j,1])-mid+n+1
    tmp_vector[pos]=as.numeric(tmp[j,2])
  }
  c_matrix[i,]=tmp_vector
}
#write.table(c_matrix, 'rpkm_c_500.txt', sep='\t', quote=F, col.names=F, row.names=F)
write.table(c_matrix, 'rpkm_200_c.txt', sep='\t', quote=F, col.names=F, row.names=F)

# load the trans ends file
out=read.table('~/Documents/LAB/data_analysis/Transcription/TransEnds_genes.txt',header=T,stringsAsFactors=F)
mid=data.frame(matrix(NA,ncol=5,nrow=nrow(out)))
for (i in 1:nrow(out))
{
  if (out[i,4]=='+')
  {
    mid[i,]=out[i,c(1,2,6,4,5)]
  }
  if (out[i,4]=='-')
  {
    mid[i,]=out[i,c(1,6,3,4,5)]
  }
}


midmatrix=as.data.frame(matrix(NA,ncol=4,nrow=1))
for (i in 1:nrow(mid))
{
  newline=c(mid[i,1],as.numeric(mid[i,2]),as.numeric(mid[i,3]),mid[i,4])
  midmatrix=rbind(midmatrix,newline)
}
midmatrix=midmatrix[2:nrow(midmatrix),]

# ORF regions
midmatrix=cbind(midmatrix, 1000/(as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])))


# STRAND==+
midmatrix_w=midmatrix[which(midmatrix[,4]=='+'),]

# STRAND==-
midmatrix_c=midmatrix[which(midmatrix[,4]=='-'),]

# normalize genes on the watson strand
red1=data1
for (i in 1:16)
{
  red1[[i]][,1]=as.numeric(red1[[i]][,1])
  red1[[i]][,2]=as.numeric(red1[[i]][,2])
}

wigmatrix_w=matrix(nrow=nrow(midmatrix_w),ncol=1000)
for (k in 1:nrow(midmatrix_w))
{
  print (k)
  ## strand
  c=as.numeric(midmatrix_w[k,1])
  start=as.numeric(midmatrix_w[k,2])
  end=as.numeric(midmatrix_w[k,3])
  tmp = red1[[c]][which(red1[[c]][,1]>=start & red1[[c]][,1]<end),]
  tmp_vector=rep(0,1000)
  for (i in 1:nrow(tmp))
  {       
    pos=floor((tmp[i,1]-start)*midmatrix_w[k,5])
    tmp_vector[pos]=as.numeric(tmp[i,2])
   }
  wigmatrix_w[k,]=tmp_vector   
   # smooth
   #bp=5
   #sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=1000/bp))
   #for (i in seq(1,1000, bp))
   #{
   #  sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
   #  sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
   #}
   #alldata_mid_w[[k]]=sm_matrix
   #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
}

# normalize genes on the crick strand
red1=data2
for (i in 1:16)
{
  red1[[i]][,1]=as.numeric(red1[[i]][,1])
  red1[[i]][,2]=as.numeric(red1[[i]][,2])
}
wigmatrix_c=matrix(nrow=nrow(midmatrix_c),ncol=1000)
for (k in 1:nrow(midmatrix_c))
{
  print (k)
  ## strand
  c=as.numeric(midmatrix_c[k,1])
  start=as.numeric(midmatrix_c[k,2])
  end=as.numeric(midmatrix_c[k,3])
  tmp = red1[[c]][which(red1[[c]][,1]>start & red1[[c]][,1]<=end),]
  tmp_vector=rep(0,1000)
  for (i in 1:nrow(tmp))
  {       
    pos=ceiling((tmp[i,1]-start-1)*midmatrix_c[k,5])
    tmp_vector[pos]=as.numeric(tmp[i,2])
  }
  wigmatrix_c[k,]=tmp_vector   
  # smooth
  #bp=5
  #sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=1000/bp))
  #for (i in seq(1,1000, bp))
  #{
  #  sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
  #  sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
  #}
  #alldata_mid_w[[k]]=sm_matrix
  #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
}



# Combine both strands
wigmatrix_c_rev=wigmatrix_c[,1000:1]

wigmatrix_w_ave = colMeans(wigmatrix_w)
wigmatrix_c_ave = colMeans(wigmatrix_c_rev)
wigmatrix_ave=(wigmatrix_w_ave+wigmatrix_c_ave)/2

png('Transcripts_ave－recom.png',width=800, height = 600)
plot(wigmatrix_ave,pch=16,col='blue',lwd=3,frame.plot=F,xlab='Normalized Transcripts',ylab='Expression',
     type='l')
dev.off()

wig_sort=cbind(1:1000, wigmatrix_ave)
# smooth
bp=5
sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=1000/bp))
for (i in seq(1,1000, bp))
{
  sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
  sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
}
png('Transcripts_ave_sm5.png',width=800, height = 600)
plot(sm_matrix,pch=16,col='blue',lwd=3,frame.plot=F,xlab='Normalized Transcripts',ylab='Expression',type='l')
dev.off()


w_ave=colMeans(w_matrix)
c_ave=colMeans(c_matrix)
wc_ave = (rev(c_ave)+w_ave)/2
png("test.png", width=800, height = 500)
plot(c(1:400)-200, wc_ave, pch=16, col='blue', frame.plot=F, xlab='Transcript Ends', ylab='Expression')
dev.off()

# normalize
#w_matrix=read.table('rpkm_w_500.txt')
w_matrix=read.table('rpkm_200_w.txt')
w_matrix_norm = w_matrix
for (i in 1:nrow(w_matrix))
{
  if (w_matrix[i,1]!=0)
  {
    x=as.vector(t(w_matrix[i,]))
    y=(x-mean(x))/sd(x)
    w_matrix_norm[i,]=y
  }
}
#write.table(w_matrix_norm,'w_matrix_500_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)
write.table(w_matrix_norm,'rpkm_200_w_norm.txt',col.names=F, row.names=F, sep='\t', quote=F)

# sort by coverage
w_matrix_norm=read.table('rpkm_200_w.txt')
#w_matrix_norm=read.table('rpkm_w_500.txt')
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(w_matrix_norm, rowSums(w_matrix_norm))
o1<- rev(order(m.row.sum[,401]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,5, length=100),seq(5,20,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:400], cluster_rows = T, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

# log rpkm
w_matrix_norm=read.table('rpkm_200_w.txt')
log_w = log2(w_matrix)
w_matrix_norm = log_w
range(w_matrix_norm)
png("test.png", width=700, height = 1200)
m.row.sum<- cbind(w_matrix_norm, rowSums(w_matrix_norm))
o1<- rev(order(m.row.sum[,401]))
m.row.sum<- m.row.sum[o1,]
bk = unique(c(seq(0,5, length=100),seq(5,10,length=100)))
hmcols<- colorRampPalette(c("white","red"))(length(bk)-1)
pheatmap(m.row.sum[,1:400], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=T, show_rownames=FALSE, show_colnames=FALSE)
dev.off()




