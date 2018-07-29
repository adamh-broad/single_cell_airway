
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

extract.field=function(string,field=1,delim="_", fixed=T) {
    return(strsplit(string,delim, fixed=fixed)[[1]][field])
}


# Logging utility function
info <- function(text, ...)
{
	cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

# Logging utility function
warn <- function(text, ...)
{
	cat(sprintf(paste(Sys.time(),"WARN:", text,"\n")))
}

### Compute TPM expression values from raw UMI counts
tpm <- function(counts, mult=10000)
{
	info("Running TPM normalisation")
	total.counts = Matrix::colSums(counts)
	scaled.counts = t(t(counts) / total.counts) 
	scaled.counts * mult
}


### Extract variable genes - custom code
get.variable.genes.umis <- function(umi.cts, residual.threshold=-0.25, use.spline=F, batch=NULL, ret.plot=F)
{
    library(Matrix)
    library(mgcv)
    if(!is.null(batch))
    {
        v = as.vector(table(batch))
        total_transcripts = data.frame(as.matrix(t(Matrix.utils::aggregate.Matrix(t(umi.cts), groupings = batch, fun="sum"))))
        detection_frac = Matrix.utils::aggregate.Matrix(t(umi.cts > 0), groupings = batch, fun="sum")
        detection_frac = data.frame(as.matrix(t(detection_frac / v)))
        test_genes = rownames(detection_frac)[Matrix::rowSums(detection_frac > 0) == length(unique(batch))]
        detection_frac = detection_frac[test_genes, ]
        total_transcripts = total_transcripts[test_genes, ]
        detection_frac$gene = rownames(detection_frac)
        total_transcripts$gene = rownames(total_transcripts)
        detection_frac = melt(detection_frac, id.vars="gene")
        colnames(detection_frac) = c("gene", "batch", "alpha")
        total_transcripts = melt(total_transcripts, id.vars="gene")
        colnames(total_transcripts) = c("gene", "batch", "UMIs")
        z = cbind(total_transcripts, detection_frac)[, c("gene", "batch", "alpha", "UMIs")]
        info("Fitting logistic GLM (controlling for batch covariate)")
        model.logit = glm(data = z, formula = alpha ~ log10(UMIs) + batch, family = binomial)
        #model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
        info("Fitting spline quantile regression (controlling for batch covariate)")
        model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15) + batch, tau=0.8)
        #model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
    }else{
        info("Computing gene dection rates (alphas)..")
        z = data.frame(UMIs = Matrix::rowSums(umi.cts), alpha= Matrix::rowSums(umi.cts>0) / ncol(umi.cts))
        z = subset(z, UMIs > 0 | alpha > 0)
        info("Fitting GLMs..")
        model.logit = glm(data = z, formula = alpha ~ log10(UMIs), family = binomial)
        #model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
        model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15), tau=0.8)
        #model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
    }
    
    
    if(use.spline){
        info("use.spline is ON. Using GAM fit (blue), logit in red")
        z$predicted = predict(object = model.gam, z, type="response")
        z$predicted.alternate = predict(object = model.logit, z, type="response")
        z$residual = model.gam$residuals
    }else{
        info("use.spline is OFF. Using logit fit (blue), GAM in red")
        z$predicted = predict(model.logit, type="response") #predict(object = model.logit, z, type="response")
        z$predicted.alternate = predict(object = model.gam, z, type="response")
        z$residual = model.logit$residuals
    }
    if(is.null(batch)) {z$gene = rownames(z)}
    outliers = subset(z, residual < residual.threshold)
    g = ggplot(z, aes(x=log10(UMIs), y=alpha, label=gene)) + geom_point(color="grey50", size=0.5, stroke=0) + 
        ylim(c(0,1)) + geom_line(aes(y=predicted), size=0.5, color="blue", linetype="dotted") + 
        geom_line(aes(y=predicted.alternate), size=0.5, color="red", linetype="dotted") + 
        geom_text(data=outliers, color="black", size=1.5, vjust=2)
    if(!is.null(batch)){g = g + facet_wrap(~batch)}
    if(!ret.plot){print(g)}
    rv = unique(unlist(lapply(rownames(outliers), extract.field, 1, delim="_")))
    if(ret.plot){return(list("var.genes"=rv, "plot"=g, "fit.data"=z, "logit"=model.logit))}
    rv
}

info("Loading default colors")
default.cols = function(n){
    #info(sprintf("Getting %s default colors", n))
    if(n<=20){
        #print(n)
        #info("Using 'Kelly' cols")
        kelly.cols(n)
    }else{
        warn("More than 20 requested, using 'Distinct' cols")
        distinct.cols(n)
    }
} 

kelly.cols <- function(n)
{
    if(n <= 20)
    {
        return(kelly[1:n])
    }else
    {
        warn("Only 20 kelly colours available")
        return(kelly)
    }
}

kelly = c(
    "#00538A", # Strong Blue
    "#C10020", # Vivid Red
    "#007D34", # Vivid Green
    "#FFB300", # Vivid Yellow
    "#803E75", # Strong Purple
    "#FF6800", # Vivid Orange
    "#A6BDD7", # Very Light Blue
    "#CEA262", # Grayish Yellow
    "#817066", # Medium Gray
    "#F6768E", # Strong Purplish Pink
    "#FF7A5C", # Strong Yellowish Pink
    "#53377A", # Strong Violet
    "#FF8E00", # Vivid Orange Yellow
    "#B32851", # Strong Purplish Red
    "#F4C800", # Vivid Greenish Yellow
    "#7F180D", # Strong Reddish Brown
    "#93AA00", # Vivid Yellowish Green
    "#593315", # Deep Yellowish Brown
    "#F13A13", # Vivid Reddish Orange
    "#232C16" # Dark Olive Green
    )

# https://www.materialui.co/colors
material.cols <- c("#f44336", #red
    "#E91E63", #pink
    "#9C27B0", #purple
    "#673AB7", #deep purple 
    "#3F51B5",  # indigo
    "#2196F3", # blue
    "#03A9F4", # light blue
    "#00BCD4", #cyan
    "#009688", # teal
    "#4CAF50", #green
    "#8BC34A", #light green
    "#CDDC39", # lime
    "#FFEB3B", #yellow
    "#FFC107", # amber
    "#FF9800", # organe
    "#FF5722", # deep orange
    "#795548", #brown
    "#9E9E9E", # grey
    "#607D8B" #blue grey
    )

# To avoid confounding of cell-type scores with cell-type quality (library complexity)
# we subtract a "control" score which is generated by averaging over a control 
# gene set. Control gene sets are chosen to contain 10 times more genes than the 
# real gene set (analogous to averaging over 10 control sets of similar size) and 
# to have the same distribution of population/bulk - based expression levels as the 
# real gene set, such that they are expected to have the same number of "zeros" and 
# to eliminate the correlation with complexity.
# ---------------------------------------------------------------------------------
get_controls <- function(counts, gene.list, verbose=F, control.genes.per.gene=10)
{
    # Find control genes by finding the closest genes in terms of expression level and % of the time we observe it
    if(verbose){info(sprintf("Finding %s background genes based on similarity to given gene set [%s genes]", 
        control.genes.per.gene*length(gene.list), length(gene.list)))}
    info("Summarizing data")
    summary = data.frame(gene=row.names(counts), mean.expr = Matrix::rowMeans(counts), fract.zero = Matrix::rowMeans(counts==0), stringsAsFactors = F)
    summary$mean.expr.s = scale(summary$mean.expr)
    summary$fract.zero.s = scale(summary$fract.zero)
    actual.genes = summary[summary$gene %in% gene.list,]
    background.genes = summary[!summary$gene %in% gene.list,]
    #find the 10 closest genes to each cell cycle marker gene and add them to the lists of control genes
    get_closest_genes <- function(i)
    {
        background.genes$dist = sqrt((background.genes$mean.expr.s - actual.genes$mean.expr.s[i])^2 + 
            (background.genes$fract.zero.s - actual.genes$fract.zero.s[i])^2)
        ordered = background.genes$gene[order(background.genes$dist)]
        ordered = ordered[!ordered %in% controls] # don't take genes that already appear in the list 
        closest = head(ordered, n=control.genes.per.gene)
        return(closest)
    }
    controls = c();
    for (i in 1:length(gene.list)){
        #info(sprintf("Finding %s control genes for %s", control.genes.per.gene, gene.list[i]))
        closest = get_closest_genes(i)
        #info(sprintf("Found %s: ", length(closest)))
        controls = unique(c(controls, closest))
    }
    if(verbose){info(sprintf("Control gene selection complete. %s genes found.", length(controls)))}
    return(controls)
}

material.heat <- function(n)
{		
		mh = c(
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        "#00BCD4", #cyan
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

### Use a generalized linear mixed model (GLMM) to 
### construct an overall estimate of a proportion, 
### (usually of lineage-labeled cells) given multiple estimates
meta_prop <- function(x, n){
    library(metafor)
    library(meta)
    if(all(x==0)){return(list("x"=NA, "ci95.high"=NA, "ci95.low"=NA))}
    if(any(n==0)){warn("Removing zero n samples (studies)"); non_zero_n = which(n>0); n = n[non_zero_n]; x = x[non_zero_n]}
    info(sprintf("Fitting.. [x=%s, n=%s]", paste(x, collapse=", "), paste(n, collapse=", ")))
    m = meta::metaprop(x, n, method = "GLMM")
    row.num = 5 #5 for random effects, 4 for fixed
    prop = as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split=" ")[[1]][8])
    ci95.low = as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split = "[", fixed=T)[[1]][2]) 
    ci95.high =  as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split = "]", fixed=T)[[2]][1])
    list("x"=prop, "ci95.high"=ci95.high, "ci95.low"=ci95.low)
}


### for each comparison, compute the y co-ordinate where the significance bar should go:
compute_signif_height <- function(an_dt, full_dt, adj = 0.02)
{
    rv = vector()
    for(i in 1:nrow(an_dt))
    {
        t1 = an_dt$start[i]
        t2 = an_dt$end[i]
        cl = an_dt$cell[i]
        rv[i] = adj + max(subset(full_dt, timepoint %in% c(t1, t2) & cell == cl)$value)
    }
    rv
}

### function to plot fit of regression
make_plot <- function(real, lin_model, quant_reg_model=NULL, title=NULL, pt.col="black", rq.se = "boot", rq.bs="xy", rq.ci=F, n_interp=100)
{

    if(!is.null(lin_model$offset)){stop("Never finished the code for a non-zero offset term!!")}
    interpolated_tps <- data.frame(timepoint = seq(min(real$timepoint), max(real$timepoint), length.out = n_interp))
    if(!is.null(lin_model$offset)) {interpolated_tps$total = n_interp}
    fitted = data.frame(predict(lin_model, interpolated_tps, type = "response", se.fit=TRUE))
    fitted <- within(fitted, {
        green_frac = fit 
        lower <- fit - 1.96 * se.fit
        higher <- fit + 1.96 * se.fit
    })
    if(!is.null(lin_model$offset)) {fitted = fitted / n_interp}
    fitted = cbind.data.frame(fitted, interpolated_tps)
    if(is.null(title)){title="Linear model"}
    g = ggplot(fitted, aes(x=timepoint, y=green_frac)) + geom_line(color="black", size=0.25, aes(y=fit), linetype="dotted") + 
        geom_point(data = real, color=pt.col) + ylab("Proportion of GFP+ cells") + xlab("") +  ggtitle(title) 
    text_label = sprintf("rate=%s", signif(lin_model$coefficients[["timepoint"]],3))
    if(!is.null(quant_reg_model))
    {
        intercept = quant_reg_model$coefficients[1]
        slope = quant_reg_model$coefficients[2]
        g = g + geom_abline(aes(slope=slope, intercept=intercept), size=0.5, color="black", alpha=.75)
        if(rq.ci){
            fitted_qr = data.frame(predict(quant_reg_model, interpolated_tps)) #, interval = "confidence", type="percentile", se="boot")) 
            fitted_qr = cbind.data.frame(fitted_qr, interpolated_tps)
            print(head(fitted_qr))
            g = g + geom_ribbon(data=fitted_qr, fill=pt.col, aes(y=coefficients, ymin=lower.bd,ymax=upper.bd),alpha=0.3, linetype="dotted") 
        }
        text_label = sprintf("%s\nrate=%s (qr)", text_label, signif(qr_lin_model$coefficients[2],3))
    }
    g = g + annotate("text", x=50, y=0.1, label= text_label, size=3)
    return(g)
}



### Estimate differences in cell-type turnover from basal cells by assessing differences 
### in the slope of the linear or quantile regression
### this pacakge (emmeans, formally lsmeans) is very useful, but doesn't work for quantreg fits
#=============================================================================================
# for the linear model:
#library(emmeans)
#m.interaction <- lm(green_frac ~ timepoint*celltype, data = d)
#m.list = emtrends(m.interaction, "celltype", var="timepoint")
#pairs(m.list)

## Emmeans doesn't work for the quantile regression model :(
#m.interaction <- quantreg::rq(formula = green_frac ~ timepoint*celltype, data = d, tau = 0.5)
#m.list = emtrends(m.interaction, "celltype", var="timepoint")
#pairs(m.list)
#=============================================================================================
## Instead, we can use the 'rank' test to compare the slopes using an interaction
## see discussion here: https://stats.stackexchange.com/questions/33013/what-test-can-i-use-to-compare-slopes-from-two-or-more-regression-models
test_pair <- function(d_all, c1, c2, test="rank", se = "boot", linear=F)
{
    ### prepare data
    x = data.frame(subset(d_all, celltype %in% c(c1, c2)))
    x$cb = paste(x$celltype, x$batch, sep="_")
    totals = aggregate(x$value, by=list(x$cb), FUN=sum)
    totals = totals[rep(1:nrow(totals),each=2),] 
    x$total = totals$x
    x = subset(x, color=="GFP")
    x$green_frac = x$value / x$total
    x$green_frac[!is.finite(x$green_frac)] <- 0
    x$celltype = factor(x$celltype)
    #print(table(x$celltype))
    if(linear)
    {
        ### fit LM models and run ANOVA LRT
        lin <- lm(green_frac ~ timepoint*celltype, data = x)
        #m.list = emmeans::emtrends(lin, "celltype", var="timepoint")
        m.list = lsmeans::lstrends(lin, "celltype", var="timepoint")
        print("regression coefficients (including slopes):")
        print(lin$coefficients)
        av = anova(lin)
        p = pairs(m.list)
        print("emtrends Comparison of Slopes:")
        print(p)
        s = summary(p)
        return(list("anova"=av, "p.value"=as.numeric(s$p.value)))
    }else{
        #print("Comparing quantreg slopes")
        ### fit quantreg models and run ANOVA LRT
        m = rq(formula = green_frac ~ timepoint, data=x, tau=0.5)
        m2 = rq(formula = green_frac ~ timepoint * celltype, data=x, tau=0.5)
        av = anova.rq(m, m2, test = test, se = se, joint=F)
        return(list("anova"=av, "p.value"=as.numeric(av$table$pvalue)))
    }
}

# utility function to generate a list of adjacent pairs of elements of a list
# c(a,b,c,d) -> list(c(a,b), c(b,c), c(c,d))
to_adjacent_pairs <- function(v){split(embed(v, 2)[, 2:1], seq(nrow(embed(v, 2)[, 2:1])))}
to_all_pairs <- function(v){y = t(combn(v, 2)); split(y, seq(nrow(y)))}


## compute and return an annotation dataframe that ggsignif can use to show p-values on the plot
## each row in the df is a p-value bar thing |-----| like this:
# #>   color start       end   y   label
# > 1     E  Good Very Good 3.6 Comp. 1
# > 2     H  Fair      Good 4.7 Comp. 2
## see here: https://www.rdocumentation.org/packages/ggsignif/versions/0.4.0
# show.pairs = list(c("Goblet", "Ciliated"), c("Ciliated", "Ionocyte"), c("Ionocyte", "Tuft"), c("Tuft", "Club"), c("Club", "Neuroendocrine")
compute_pvals <- function(d, show.pairs, se.use = "boot", test.use="rank", linear=F)
{
    rv = data.frame(t(data.frame(show.pairs, stringsAsFactors = F)), stringsAsFactors = F)
    rownames(rv) = 1:nrow(rv)
    colnames(rv) = c("start", "end")
    for(i in 1:length(show.pairs))
    {
        p = show.pairs[[i]]
        rv$p.value[i] <- test_pair(d, p[1], p[2], test=test.use, linear=linear, se=se.use)$p.value
    }
    rv$p.value = as.numeric(rv$p.value)
    rv
}


### for each comparison, compute the y co-ordinate where the significance bar should go:
compute_signif_height_slope <- function(an_dt, full_dt, adj = 0.02)
{
    rv = vector()
    for(i in 1:nrow(an_dt))
    {
        c1 = an_dt$start[i]
        c2 = an_dt$end[i]
        c1_index = which(full_dt$celltype == c1)
        c2_index = which(full_dt$celltype == c2)
        rv[i] = adj + max(full_dt$coeff_upr[c1_index], full_dt$coeff_upr[c2_index])
    }
    rv
}




