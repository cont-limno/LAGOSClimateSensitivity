# single variable plot of regions side by side
plot_var = 'summer_ppt'
gradient_var = 'wyppt_PRISM'
par(mfrow=c(1,2))
plot(PRISM_gradient_df_UM[,plot_var]~PRISM_gradient_df_UM[,gradient_var], 
     main=paste0('Midwest'), pch=20, col='red',
     xlim=c(600, 1600), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_UM[,plot_var]~ PRISM_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(PRISM_gradient_df_NE[,plot_var]~PRISM_gradient_df_NE[,gradient_var], 
     main=paste0('Northeast'), pch=20, col='dodgerblue',
     xlim=c(600,1600), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_NE[,plot_var]~ PRISM_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)
cor(PRISM_gradient_df$wyppt_PRISM, PRISM_gradient_df$summer_ppt)