# Import libraries
library(sf)
library(tidyverse)
library(GWmodel)    # to undertake the GWR
library(tmap)       # for mapping
library(spdep)
library(ggplot2)
library(caret)
library(ggpubr)
library(gwrr)

# Import Map
Jateng_map <- read_sf('RBI_50K_2023_Jawa Tengah.x26272/RBI_50K_2023_Jawa Tengah.shp')
names(Jateng_map)
Jateng_map <- Jateng_map[,c(1,26)]
colnames(Jateng_map) <- c("Kab","geometry")
#plot(Jateng_map)
Jateng_map$Kab

# Import Data
library(readxl)
df<- read_excel('Koordinat Jawa Tengah.xlsx')
names(df)
df <- df[,-c(2:3)]
head(df)

# Analisis Deskriptif
summary_data <- df[,-1]
summary_df <- data.frame(
  Variable = names(summary_data),
  Min = apply(summary_data, 2, min),
  Median = apply(summary_data, 2, median),
  Mean = apply(summary_data, 2, mean),
  Max = apply(summary_data, 2, max)
)

#writexl::write_xlsx(summary_df,"Analisis Deskriptif.xlsx") 

# Scaling
df[,-1] <- scale(df[,-1])
Data <- merge(Jateng_map,df,by="Kab") # Gabungkan data peta dan data tabular/excel
Data <- st_as_sf(Data)
Data <- st_zm(Data)

ggplot(Data) + geom_sf(aes(fill = Y)) +
  scale_fill_gradient2(
    midpoint = mean(Data$Y), low = "#28E2E5", mid = "white", high = "#DF536B"
  ) +
  theme_bw()+
  geom_text(
    aes(label = Kab, x = coordinates(as(Data,"Spatial"))[,1], y = coordinates(as(Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = FALSE
  )+ggtitle("Sebaran Data Y di Jawa Tengah")+xlab("Longitude")+ylab("Latitude")

# OLS Model
formula <- paste("Y ~ X2+X3+X4+X5+X6+X7+X8+X9+X10")
ols = lm(formula, data = df[,-1])
ols_summary <- summary(ols)
ols_summary

# determine studentised residuals and attach to data
s.resids <- rstudent(ols)

# Normality test
shapiro.test(s.resids)

# Heteroskedastisity test
lmtest::bptest(ols)

# Multicolinierity test
library(car)
vif(ols)

# Moran test
nb <- poly2nb(st_zm(Jateng_map), queen=TRUE)
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
moran.test(Data$Y,lw, alternative="greater")

# GWR
# convert to sp
Data.sp = as(Data, "Spatial")
# determine the kernel bandwidth
bw <- bw.gwr(formula,
             approach = "AIC",
             adaptive = T,
             data=Data.sp)
bw

# fit the GWR model
gwr.model <- gwr.basic(formula,
                            adaptive = T,
                            data = Data.sp,
                            kernel = "bisquare",
                            bw = bw)
gwr_local_vif<- gwr.collin.diagno(formula,
                  adaptive = T,
                  data = Data.sp,
                  kernel = "bisquare",
                  bw = bw)
gwr_local_vif$VIF
gwr_sf = st_as_sf(gwr.model$SDF)

# Menampilkan ringkasan model
summary(gwr_sf)

# GWRR
# fit the GWR model
locs <- cbind(coordinates(Data.sp)[,1],coordinates(Data.sp)[,2])
gwrr.model <- gwrr.est(Y ~ X2+X3+X4+X5+X6+X7+X8+X9+X10, locs, data.frame(Data.sp), "exp",bw=TRUE,rd=TRUE,cv.tol = 0.000001)

gwrr_local_vif<- gwr.collin.diagno(formula,
                  adaptive = T,
                  data = Data.sp,
                  kernel = "exponential",
                  bw = gwrr.model$phi)
gwrr_local_vif$VIF

# Menampilkan ringkasan model
gwrr.model

# Menghitung nilai VIF dari model GWRR
# Loop untuk setiap lokasi i
X <- as.matrix(df[,3:12])
n <- nrow(X)
H_local <- matrix(0, nrow = n, ncol = ncol(X))
for (i in 1:n) {
  Wi <- diag(W_GWRR[i,]) # Matriks bobot spasial untuk lokasi i
  Hi <- X %*% solve(t(X) %*% Wi %*% X) %*% t(X) %*% Wi # Matriks hat lokal untuk lokasi i
  H_local[,i] <- diag(Hi) # Elemen diagonal
}

# Hitung VIF lokal
VIF_local <- 1 / (1 - H_local)

# Tampilkan VIF lokal (perlu pengolahan lebih lanjut untuk presentasi yang baik)
print(VIF_local)

# Pemilihan Model Terbaik
rmse <- function(residual){
  sqrt(mean((residual)^2))
}

overall_summary <- data.frame(R2 = c(ols_summary$r.squared,gwr.model$GW.diagnostic$gw.R2,gwrr.model$rsquare),
                              RMSE = c(rmse(ols_summary$residuals),rmse(gwr_sf$residual),gwrr.model$RMSE))
rownames(overall_summary) <- c("OLS","GWR","GWRR")
overall_summary

writexl::write_xlsx(overall_summary,"Perbandingan antar Model.xlsx")

## Menghitung signifikansi tiap variabel
# Menghitung standar deviasi dari koefisien lokal
coef_sd <- apply(t(gwrr.model$beta), 2, sd)

# Menghitung interval kepercayaan 95%
alpha <- 0.05
z_critical <- qnorm(1 - alpha / 2)  # Nilai kritis untuk 95%

ci_lower <- t(gwrr.model$beta) - z_critical * coef_sd
ci_upper <- t(gwrr.model$beta) + z_critical * coef_sd

# Cek apakah nol berada dalam interval kepercayaan
significant <- (ci_lower > 0) | (ci_upper < 0)

# Menampilkan variabel yang signifikan secara spasial
significant <- ifelse(significant==TRUE,"Significant","Not Significant")
colnames(significant) <- c("signif_Intercept","signif_X1","signif_X2","signif_X3","signif_X4","signif_X5","signif_X10","signif_X15")
significant <- data.frame(significant)
significant$Kab <- Data$Kab
significant

All_Data <- merge(st_zm(Data),significant,by = "Kab")
writexl::write_xlsx(All_Data,"Signifikansi Tiap Variabel.xlsx")

# Plot Signifikansi Hasil Model Terbaik
# X1
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X1)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X1")+xlab("Longitude")+ylab("Latitude")

# X2
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X2)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X2")+xlab("Longitude")+ylab("Latitude")

# X3
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X3)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X3")+xlab("Longitude")+ylab("Latitude")

# X4
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X4)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X4")+xlab("Longitude")+ylab("Latitude")

# X5
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X5)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X5")+xlab("Longitude")+ylab("Latitude")

# X10
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X10)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X10")+xlab("Longitude")+ylab("Latitude")

# X15
ggplot(data=All_Data) +
  geom_sf(mapping=aes(fill =signif_X15)) +
  scale_fill_manual(values = c("#DF536B","#28E2E5"))+
  labs(fill="Significancy")+
  geom_text(
    aes(label = Kab, x = coordinates(as(All_Data,"Spatial"))[,1], y = coordinates(as(All_Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X15")+xlab("Longitude")+ylab("Latitude")
