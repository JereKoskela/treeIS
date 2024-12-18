library("viridis")

# =================== cost with 100k particles under HUW =======================
cost_55 = read.table("./cost-55-huw.txt")
cost_550 = read.table("./cost-550-huw.txt")
cost_5500 = read.table("./cost-5500-huw.txt")
ymax = log10(max(cost_55, cost_550, cost_5500))
col = viridis(3)
pdf("cost-100k-huw.pdf")
plot(cost_55[,1], log10(cost_55[,2]), type="l", lwd=2, xlab = "Fraction of sample size", ylab = "log(E[W^2] / E[W]^2) (base 10)", cex.lab=1.3, cex.axis=1.3, ylim = c(0, ymax), main = "Hobolth-Uyenoyama-Wiuf", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(cost_550[,1], log10(cost_550[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[2])
par(new=TRUE)
plot(cost_5500[,1], log10(cost_5500[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[3])
legend("topright", legend=c("n = 55 (6s)", "n = 550 (1m 35s)", "n = 5500 (18m 49s)"), lwd=c(2,2,2), col = col, cex=1.3)
dev.off()
# ==============================================================================

# =================== cost with 100k particles under SD ========================
cost_55 = read.table("./cost-55-sd.txt")
cost_550 = read.table("./cost-550-sd.txt")
cost_5500 = read.table("./cost-5500-sd.txt")
col = viridis(3)
ymax = log10(max(cost_55, cost_550, cost_5500))
pdf("cost-100k-sd.pdf")
plot(cost_55[,1], log10(cost_55[,2]), type="l", lwd=2, xlab = "Fraction of sample size", ylab = "log(E[W^2] / E[W]^2) (base 10)", cex.lab=1.3, cex.axis=1.3, ylim = c(0, ymax), main = "Stephens-Donnelly", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(cost_550[,1], log10(cost_550[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[2])
par(new=TRUE)
plot(cost_5500[,1], log10(cost_5500[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[3])
legend("topright", legend=c("n = 55 (2s)", "n = 550 (31s)", "n = 5500 (6m 14s)"), lwd=c(2,2,2), col = col, cex=1.3)
dev.off()
# ==============================================================================

# =================== cost with 100k particles under GT ========================
cost_55 = read.table("./cost-55-gt.txt")
cost_550 = read.table("./cost-550-gt.txt")
cost_5500 = read.table("./cost-5500-gt.txt")
col = viridis(3)
ymax = log10(max(cost_55, cost_550, cost_5500))
pdf("cost-100k-gt.pdf")
plot(cost_55[,1], log10(cost_55[,2]), type="l", lwd=2, xlab = "Fraction of sample size", ylab = "log(E[W^2] / E[W]^2) (base 10)", cex.lab=1.3, cex.axis=1.3, ylim = c(0, ymax), main = "Griffiths-Tavare", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(cost_550[,1], log10(cost_550[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[2])
par(new=TRUE)
plot(cost_5500[,1], log10(cost_5500[,2]), type="l", lwd=2, xlab = "", ylab = "", xaxt="n", yaxt="n", ylim = c(0, ymax), col=col[3])
legend("topright", legend=c("n = 55 (6s)", "n = 550 (1m)", "n = 5500 (22m 7s)"), lwd=c(2,2,2), col = col, cex=1.3)
dev.off()
# ==============================================================================


# ========= Data for infinite sites likelihood surfaces ========================
likelihood_100k_55 = read.table("./likelihood-100k-results-55-huw.txt")
likelihood_200k_550 = read.table("./likelihood-200k-results-550-huw.txt")
likelihood_variable_55 = read.table("./likelihood-variable-results-55-huw.txt")
likelihood_variable_550 = read.table("./likelihood-variable-results-550-huw.txt")
likelihood_1000_55 = read.table("./likelihood-1000-results-55-huw.txt")
likelihood_2000_550 = read.table("./likelihood-2000-results-550-huw.txt")
likelihood_mean_55 = read.table("./likelihood-mean-results-55-huw.txt")
likelihood_mean_550 = read.table("./likelihood-mean-results-550-huw.txt")

# ======== n = 55 likelihood surface as per Hobolth et al (2008) ===============
th = seq(2, 10, by=1)
p = 20
col=viridis(4)
sch4_sample = c(37000, 49600, 58600, 64000, 69400, 73000, 74800, 76600, 78400)
pdf("likelihood-55-huw.pdf")
plot(th, exp(likelihood_100k_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 55 (HUW)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (1m 16s)", "Schedule 2 (54s)", "Schedule 3 (1s)", "Schedule 4 (50s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ======== n = 550 likelihood surface as per Hobolth et al (2008) ==============
th = seq(2, 10, by=1)
p = 46
col=viridis(4)
sch4_sample = c(71120, 97760, 116120, 128720, 138440, 145640, 151400, 156080, 159680)
pdf("likelihood-550-huw.pdf")
plot(th, exp(likelihood_200k_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 550 (HUW)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (1h 17m)", "Schedule 2 (46m 49s)", "Schedule 3 (39s)", "Schedule 4 (52m 25s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ====== Data for HUW infinite sites likelihood surfaces with resampling========
likelihood_100k_55 = read.table("./likelihood-100k-results-55-huw-re.txt")
likelihood_200k_550 = read.table("./likelihood-200k-results-550-huw-re.txt")
likelihood_variable_55 = read.table("./likelihood-variable-results-55-huw-re.txt")
likelihood_variable_550 = read.table("./likelihood-variable-results-550-huw-re.txt")
likelihood_1000_55 = read.table("./likelihood-1000-results-55-huw-re.txt")
likelihood_2000_550 = read.table("./likelihood-2000-results-550-huw-re.txt")
likelihood_mean_55 = read.table("./likelihood-mean-results-55-huw-re.txt")
likelihood_mean_550 = read.table("./likelihood-mean-results-550-huw-re.txt")

# ======== n = 55 likelihood surface as per Hobolth et al (2008) ===============
th = seq(2, 10, by=1)
p = 20
col=viridis(4)
sch4_sample = c(37000, 49600, 58600, 64000, 69400, 73000, 74800, 76600, 78400)
pdf("likelihood-55-huw-re.pdf")
plot(th, exp(likelihood_100k_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 55 (HUW + Re)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (1m 34s)", "Schedule 2 (47s)", "Schedule 3 (1s)", "Schedule 4 (1m 8s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ======== n = 550 likelihood surface as per Hobolth et al (2008) ==============
th = seq(2, 10, by=1)
p = 46
col=viridis(4)
sch4_sample = c(71120, 97760, 116120, 128720, 138440, 145640, 151400, 156080, 159680)
pdf("likelihood-550-huw-re.pdf")
plot(th, exp(likelihood_200k_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 550 (HUW + Re)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (1h 15m)", "Schedule 2 (46m 12s)", "Schedule 3 (45s)", "Schedule 4 (42m 22s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ========= Data for SD infinite sites likelihood surfaces =====================
likelihood_100k_55 = read.table("./likelihood-100k-results-55-sd.txt")
likelihood_200k_550 = read.table("./likelihood-200k-results-550-sd.txt")
likelihood_variable_55 = read.table("./likelihood-variable-results-55-sd.txt")
likelihood_variable_550 = read.table("./likelihood-variable-results-550-sd.txt")
likelihood_1000_55 = read.table("./likelihood-1000-results-55-sd.txt")
likelihood_2000_550 = read.table("./likelihood-2000-results-550-sd.txt")
likelihood_mean_55 = read.table("./likelihood-mean-results-55-sd.txt")
likelihood_mean_550 = read.table("./likelihood-mean-results-550-sd.txt")

# ======== n = 55 likelihood surface as per Stephens & Donnelly (2000) =========
th = seq(2, 10, by=1)
p = 20
col=viridis(4)
sch4_sample = c(37000, 49600, 58600, 64000, 69400, 73000, 74800, 76600, 78400)
pdf("likelihood-55-sd.pdf")
plot(th, exp(likelihood_100k_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 55 (SD)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (31s)", "Schedule 2 (23s)", "Schedule 3 (1s)", "Schedule 4 (22s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ======= n = 550 likelihood surface as per Stephens & Donnelly (2000) =========
th = seq(2, 10, by=1)
p = 46
col=viridis(4)
sch4_sample = c(71120, 97760, 116120, 128720, 138440, 145640, 151400, 156080, 159680)
pdf("likelihood-550-sd.pdf")
plot(th, exp(likelihood_200k_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 550 (SD)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (12m 20s)", "Schedule 2 (9m 19s)", "Schedule 3 (9s)", "Schedule 4 (8m 20s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ====== Data for infinite sites likelihood surfaces with resampling============
likelihood_100k_55 = read.table("./likelihood-100k-results-55-sd-re.txt")
likelihood_200k_550 = read.table("./likelihood-200k-results-550-sd-re.txt")
likelihood_variable_55 = read.table("./likelihood-variable-results-55-sd-re.txt")
likelihood_variable_550 = read.table("./likelihood-variable-results-550-sd-re.txt")
likelihood_1000_55 = read.table("./likelihood-1000-results-55-sd-re.txt")
likelihood_2000_550 = read.table("./likelihood-2000-results-550-sd-re.txt")
likelihood_mean_55 = read.table("./likelihood-mean-results-55-sd-re.txt")
likelihood_mean_550 = read.table("./likelihood-mean-results-550-sd-re.txt")

# ======== n = 550 likelihood surface as per Stephens & Donnelly (2000) ========
th = seq(2, 10, by=1)
p = 20
col=viridis(4)
sch4_sample = c(37000, 49600, 58600, 64000, 69400, 73000, 74800, 76600, 78400)
pdf("likelihood-55-sd-re.pdf")
plot(th, exp(likelihood_100k_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 55 (SD + Re)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_100k_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_100k_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_55[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_1000_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_1000_55[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 15), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_55[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_55[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,15), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (41s)", "Schedule 2 (25s)", "Schedule 3 (1s)", "Schedule 4 (20s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================

# ======= n = 550 likelihood surface as per Stephens & Donnelly (2000) =========
th = seq(2, 10, by=1)
p = 46
col=viridis(4)
sch4_sample = c(71120, 97760, 116120, 128720, 138440, 145640, 151400, 156080, 159680)
pdf("likelihood-550-sd-re.pdf")
plot(th, exp(likelihood_200k_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "theta", ylab = paste0("Likelihood * 10^", p), cex.lab=1.3, cex.axis=1.3, main = "n = 550 (SD + Re)", cex.main = 1.3, col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_200k_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_200k_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[1])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_variable_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_variable_550[,3] + 2 * p * log(10) - 5 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[2])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_2000_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_2000_550[,3] + 2 * p * log(10) - 3 * log(10))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[3])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)), type = "l", lwd=2, ylim=c(0, 12), xlab = "", ylab = "", xaxt="n", yaxt="n", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) + 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
par(new=TRUE)
plot(th, exp(likelihood_mean_550[,2] + p * log(10)) - 2 * sqrt(exp(likelihood_mean_550[,3] + 2 * p * log(10) - log(sch4_sample))), type = "l", lwd=2, lty=2, ylim=c(0,12), xaxt="n", yaxt="n", xlab="", ylab="", col=col[4])
legend("topright", legend = c("Schedule 1 (10m 48s)", "Schedule 2 (8m 35s)", "Schedule 3 (6s)", "Schedule 4 (6m 31s)", "2 Standard errors"), lty = c(1, 1, 1, 1, 2), lwd=c(2,2,2,2,2), col=c(col, "black"), cex=1.3)
dev.off()
# ==============================================================================