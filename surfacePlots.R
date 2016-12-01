testing rsm for fitting and displaying surface plot
require(rsm)

https://cran.r-project.org/web/packages/rsm/vignettes/rsm-plots.pdf
https://cran.r-project.org/web/packages/rsm/vignettes/rsm.pdf


CR1 <- combined.resp[ treatment %in% "mercaptopurine" & fingerprints %in% "GFP_pos.2m_Srxn1", ]
CR1 <- CR1[, mean(value), by = c("timeID","dose_uM")]

CR1 <- CR1[ timeID < 17]
setnames(CR1, "V1", "value")

CR1.lmP <- lm(value ~ poly(dose_uM * timeID, degree= 3), data = CR1)
persp(CR1.lmP,  dose_uM ~timeID, zlab = "response", zlim= c(0,1))
anova(CR1.lmP)


CR1.lmP.dmso <- lm(value ~ poly(dose_uM * timeID, degree= 3), data = CR1)
persp(CR1.lmP.dmso,  dose_uM ~timeID, zlab = "response", zlim= c(0,1))

anova(CR1.lmP, CR1.lmP.dmso)




text()

#
col.matrix <- matrix(runif(80), nrow = 16, ncol = 5)

#CR1 <- CR1[, list( timeID, dose_uM, value )]



pdf("test.pdf", height = 60, width = 60)
par(mfrow = c(20,20))


text()
dev.off()

# conclusie:

# model per dose de time courses
# maak een grid: per compound-replicate een matrix
# maak een matrix met t-test waarden tussen de replicates
# bereken gemiddelde van de replicate gemiddelde responsen
# plot met 'persp' en kleur met de p-waarden


persp(seq(10, 300, 5), seq(10, 300, 5), z, phi = 45, theta = 45,
      xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
      main = "Surface elevation data"
)

