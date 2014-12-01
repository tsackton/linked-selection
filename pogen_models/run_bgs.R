#### define functions ###
require(data.table)
require(minpack.lm)

get_neut = function(pi) {
  mod=try(nlsLM(pi ~ neutral_pi, start=c(neutral_pi=mean(pi)), lower=c(0), upper=c(1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod),ll=logLik(mod),pi.neut=coef(mod)[1],param2=0)
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA),param2=as.numeric(NA))
  }
}

get_hh_nofd = function(pi,rr) {
  mod=try(nlsLM(pi ~ (neutral_pi*rr)/(rr + (alpha)), start=c(neutral_pi=mean(pi),alpha=1e-9), lower=c(0,0), upper=c(1,.1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod),ll=logLik(mod),pi.neut=coef(mod)[1],param2=coef(mod)[2])
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

get_bgs_nofd = function(pi,rr) {
  mod=try(nlsLM(pi ~ (neutral_pi * exp(-1*mu/rr)), start=c(neutral_pi=mean(pi),mu=1e-9), lower=c(0,0), upper=c(1,.1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod),ll=logLik(mod),pi.neut=coef(mod)[1],param2=coef(mod)[2])
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

get_bgs_fd = function(pi,gk) {
  mod=try(nlsLM(pi ~ neutral_pi*exp(-1*gk), start=c(neutral_pi=mean(pi)), lower=c(0), upper=c(1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod), ll=logLik(mod),pi.neut=coef(mod)[1],param2=0)
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

get_hh_fd = function(pi,rr,fd) {
  mod=try(nlsLM(pi ~ (neutral_pi*rr)/(rr + (alpha*fd)), start=c(neutral_pi=mean(pi),alpha=1e-9), lower=c(0,0), upper=c(1,.1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod), ll=logLik(mod),pi.neut=coef(mod)[1],param2=coef(mod)[2])
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

get_full_nofd = function(pi,gk,rr) {
  mod=try(nlsLM(pi ~ (neutral_pi*exp(-1*gk)*rr)/(rr + (alpha*exp(-1*gk))), start=c(neutral_pi=mean(pi),alpha=1e-9), lower=c(0,0), upper=c(1,.1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod), ll=logLik(mod),pi.neut=coef(mod)[1],param2=coef(mod)[2])
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

get_full_fd = function(pi,gk,rr,fd) {
  mod=try(nlsLM(pi ~ (neutral_pi*exp(-1*gk)*rr)/(rr + (alpha*fd*exp(-1*gk))), start=c(neutral_pi=mean(pi),alpha=1e-9), lower=c(0,0), upper=c(1,.1), control=nls.lm.control(maxiter=200)))
  if (!inherits(mod, "try-error")) {
    list(aic=AIC(mod), ll=logLik(mod),pi.neut=coef(mod)[1],param2=coef(mod)[2])
  } else {
    list(aic=as.numeric(NA),ll=as.numeric(NA), pi.neut=as.numeric(NA), param2=as.numeric(NA))
  }
}

recrates = fread("windows_for_theta.out", header=T);
#convert recombination rate from cM/Mb to recombination events per site per generation
recrates$rho = (1-exp(-1*(recrates$rr/1e6)*2/100))/2
recrates = subset(recrates, recrates$species == spec & recrates$ratemtd==recmethod, select=c("chr", "window.id", "rho", "pi", "fd", "use", "filt"));
recrates.use = subset(recrates,use==1)

### actual code ##

bgs<-fread(file, header=F)
setnames(bgs, c("sp", "chr", "window.id", "U", "sh", "P", "G"))
bgs.use<-merge(bgs, recrates.use, by=c("chr", "window.id"),allow.cartesian=T)
setkeyv(bgs.use, c("sh", "P", "filt"))
bgs.use=bgs.use[complete.cases(bgs.use),]

sel.full.fd<-bgs.use[,get_full_fd(pi,G,rho,fd),by="sh,P,filt"]
sel.full.fd$model="full(hh_var)"
sel.full.nofd<-bgs.use[,get_full_nofd(pi,G,rho),by="sh,P,filt"]
sel.full.nofd$model="full(hh_const)"
sel.hh.fd<-bgs.use[,get_hh_fd(pi,rho,fd),by="sh,P,filt"]
sel.hh.fd$model="hh_only(var)"
sel.hh.nofd<-bgs.use[,get_hh_nofd(pi,rho),by="sh,P,filt"]
sel.hh.nofd$model="hh_only(const)"
sel.bgs.fd<-bgs.use[,get_bgs_fd(pi,G),by="sh,P,filt"]
sel.bgs.fd$model="bgs_only(var)"
sel.bgs.nofd<-bgs.use[,get_bgs_nofd(pi,rho),by="sh,P,filt"]
sel.bgs.nofd$model="bgs_only(const)"
sel.neut<-bgs.use[,get_neut(pi),by="sh,P,filt"]
sel.neut$model="neutral"

#write out all results
outdir="modres/"
all.res<-rbind(sel.full.fd,sel.full.nofd,sel.hh.fd,sel.hh.nofd,sel.bgs.fd,sel.bgs.nofd,sel.neut)
res.file=paste(spec, as.character(windsize), as.character(uversion), as.character(recmethod), "select.out", sep="_")
res.file=paste0(as.character(outdir), res.file)
write.table(all.res, file=res.file, sep="\t", quote=F, row.names=F, col.names=T)
