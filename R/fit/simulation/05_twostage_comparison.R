## Set FIT <- FALSE to skip the (long) saemix runs and only rebuild the figure from an
## existing comparacion_twostage_vs_joint.csv.
if (!exists("FIT")) FIT <- TRUE

suppressMessages(library(saemix))
library(ggplot2)
setwd(file.path(REPO, "monolix", "data1", "Model_1"))   # REPO = path to your clone (define it once; see README)
if (FIT) {
real <- read.csv("jointmodel_data_new.csv")
ids <- sort(unique(real$id))
design <- lapply(ids, function(i){di<-real[real$id==i,]; list(tl=di$time[di$ytype==1], tb=di$time[di$ytype==2])})
truth <- c(a_pop=4.54,b_pop=15.6,c_pop=7.22,alpha_pop=34.3,beta_pop=-7.94,omega_a=0.502,omega_b=3.31,omega_c=1.95)
sr <- 0.261
sim1 <- function(seed){ set.seed(seed); n<-length(design)
  a<-rnorm(n,truth["a_pop"],truth["omega_a"]); b<-rnorm(n,truth["b_pop"],truth["omega_b"]); cc<-pmax(rnorm(n,truth["c_pop"],truth["omega_c"]),0.1)
  long<-NULL; D<-numeric(n)
  for(k in 1:n){ tl<-design[[k]]$tl; E<-a[k]/(1+exp(-(tl-b[k])/cc[k])); y<-E+rnorm(length(tl),0,sr)
    long<-rbind(long,data.frame(id=k,time=tl,y=y)); D[k]<-rbinom(1,1,1/(1+exp(-(truth["alpha_pop"]+truth["beta_pop"]*a[k])))) }
  list(long=long, D=D, a_true=a) }
modlog <- function(psi,id,xidep){ t<-xidep[,1]; a<-psi[id,1];b<-psi[id,2];c<-psi[id,3]; a/(1+exp(-(t-b)/c)) }
twostage <- function(dat){
  sd<-saemixData(name.data=dat$long,name.group="id",name.predictors="time",name.response="y",verbose=FALSE)
  sm<-saemixModel(model=modlog,psi0=matrix(c(4.54,15.6,7.22),nrow=1,dimnames=list(NULL,c("a","b","c"))),
                  transform.par=c(0,0,0),covariance.model=diag(3),error.model="constant",verbose=FALSE)
  fit<-saemix(sm,sd,list(save=FALSE,save.graphs=FALSE,print=FALSE,displayProgress=FALSE,nbiter.saemix=c(200,100)))
  ahat<-psi(fit)[,1]; g<-glm(dat$D~ahat,family=binomial)
  ph<-predict(g,type="response")
  auc<-{o<-order(ph); r<-rank(ph); n1<-sum(dat$D==1);n0<-sum(dat$D==0); (sum(r[dat$D==1])-n1*(n1+1)/2)/(n1*n0)}
  c(eta1=unname(coef(g)[2]), astar=unname(-coef(g)[1]/coef(g)[2]), auc=auc)
}
joint <- read.csv("recuperacion_parametros_M1_resultados.csv")
R<-50; out<-data.frame()
for(r in 1:R){ dat<-sim1(10000+r); ts<-tryCatch(twostage(dat),error=function(e)c(NA,NA,NA))
  out<-rbind(out,data.frame(rep=r, ts_eta1=ts[1], ts_astar=ts[2], ts_auc=ts[3],
    joint_eta1=joint$est_beta_pop[joint$rep==r],
    joint_astar=-joint$est_alpha_pop[joint$rep==r]/joint$est_beta_pop[joint$rep==r])) }
write.csv(out,"comparacion_twostage_vs_joint.csv",row.names=FALSE)
}  # end if (FIT)
out <- read.csv("comparacion_twostage_vs_joint.csv")
jdiv <- abs(out$joint_eta1) > 200
cat("\n===== COMPARACION two-stage vs joint (R=50, verdad eta1=-7.94, a*=4.32) =====\n")
cat(sprintf("two-stage eta1:  media %.2f  sd %.2f   [atenuacion: |eta1|/7.94 = %.0f%%]\n",
    mean(out$ts_eta1), sd(out$ts_eta1), 100*mean(abs(out$ts_eta1))/7.94))
cat(sprintf("joint eta1 (no-div): media %.2f  sd %.2f   (divergen %.0f%% de las reps)\n",
    mean(out$joint_eta1[!jdiv]), sd(out$joint_eta1[!jdiv]), 100*mean(jdiv)))
cat(sprintf("two-stage a*:  media %.3f  sd %.3f\n", mean(out$ts_astar), sd(out$ts_astar)))
cat(sprintf("joint a*:      media %.3f  sd %.3f\n", mean(out$joint_astar,na.rm=T), sd(out$joint_astar,na.rm=T)))
cat(sprintf("two-stage AUC (in-sample): media %.3f\n", mean(out$ts_auc)))

## --- Figure (Section 3.5): attenuation of eta_1 and agreement on a* ---
d <- out
jdiv2 <- abs(d$joint_eta1) > 200
e1 <- rbind(
  data.frame(quantity="Association slope (eta_1)", method="Two-stage", value=d$ts_eta1),
  data.frame(quantity="Association slope (eta_1)", method="Joint",     value=d$joint_eta1[!jdiv2]))
as_ <- rbind(
  data.frame(quantity="Decision boundary (a*)", method="Two-stage", value=d$ts_astar),
  data.frame(quantity="Decision boundary (a*)", method="Joint",     value=d$joint_astar))
long <- rbind(e1, as_); long$method <- factor(long$method, levels=c("Two-stage","Joint"))
truths <- data.frame(quantity=c("Association slope (eta_1)","Decision boundary (a*)"), t=c(-7.94, 4.32))
p <- ggplot(long, aes(method, value, fill=method)) +
  geom_boxplot(width=0.55, outlier.size=0.8, alpha=0.75) +
  geom_jitter(width=0.12, size=0.7, alpha=0.4) +
  geom_hline(data=truths, aes(yintercept=t), linetype="dashed", colour="red") +
  facet_wrap(~quantity, scales="free_y") +
  scale_fill_manual(values=c("Two-stage"="#66c2a5","Joint"="#fc8d62"), guide="none") +
  labs(x=NULL, y="Estimate", caption="Red dashed line: true value") +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        plot.background=element_rect(fill="white",colour=NA))
out_fig <- file.path(REPO, "figures"); dir.create(out_fig, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(out_fig,"comparison_twostage.pdf"), p, width=7.2, height=3.8)
cat("Figure written: figures/comparison_twostage.pdf\n")
