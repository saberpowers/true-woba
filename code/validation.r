require(glmnet)
require(lme4)
require(Matrix)
require(xtable)
require(weights)
source('code/util.r')




# DATA PROCESSING

# read data
data = read.csv('data/all2015.csv', header = FALSE)
colnames(data) = read.csv('data/fields.csv')$Header

# create y (vector containing outcome variable)
events = sort(c('F', 'G', 'K', 'BB', 'HBP', '1B', '2B', '3B', 'HR'))
y = as.character(data$EVENT_CD)
y[(y == 2 | y == 18 | y == 19) & data$BATTEDBALL_CD != 'G'] = 'F'
y[(y == 2 | y == 18 | y == 19) & data$BATTEDBALL_CD == 'G' & !data$SH_FL] = 'G'
y[y == 3] = 'K'
y[y == 14] = 'BB'
y[y == 16] = 'HBP'
y[y == 20] = '1B'
y[y == 21] = '2B'
y[y == 22] = '3B'
y[y == 23] = 'HR'
subset = is.element(y, events) & data$BAT_FLD_CD != 1
y = as.factor(y[subset])

# create x (_sparse_ design matrix containing covariates for regression)
hand = (data$RESP_BAT_HAND_CD != data$RESP_PIT_HAND_CD)[subset]
home = data$BAT_HOME_ID[subset] == 1
stadium = as.numeric(as.factor(as.character(data$HOME_TEAM_ID[subset]))) + 2
batter = as.numeric(as.factor(as.character(data$BAT_ID[subset]))) + max(stadium)
pitcher = as.numeric(as.factor(as.character(data$PIT_ID[subset]))) + max(batter)
x = sparseMatrix(c(which(hand), which(home), rep(1:sum(subset), 3)),
    c(rep(1, sum(hand)), rep(2, sum(home)), stadium, batter, pitcher))




# CODE TO PRODUCE FIGURE AND TABLES IN METHODS SECTION OF PAPER

# hold out test set (randomly selection half of the data)
set.seed(4987)
test = sample(1:sum(subset), floor(sum(subset)/2))

# fit naive estimator for each outcome on training data
count = aggregate(y[-test], by = list(batter[-test]), table)
pa = rowSums(count[,-1])
naive = count[,-1]/rowSums(count[,-1])
names(pa) = rownames(naive) = sort(unique(data$BAT_ID[subset][-test]))

# fit null estimator for each outcome on training data
null = naive
for (e in events) {
    null[, e] = mean(y[-test] == e)
}

# fit regressed estimator for each outcome on training data
regressed = naive
for (e in events) {
    priorMean = mean(y[-test] == e)
    var = priorMean*(1-priorMean)/pa
    priorVar = estimatePopulationVariance(naive[, e], var, priorMean)
    regressed[, e] = (naive[, e]/var + priorMean/priorVar)/(1/var + 1/priorVar)
}
regressed = regressed/rowSums(regressed)

# fit ridge regression estimator for each outcome on training data
Sys.time = Sys.time()
batter.b = as.numeric(as.factor(as.character(data$BAT_ID[subset][-test])))
x.b = sparseMatrix(1:length(batter.b), batter.b)
path = exp(seq(-6, -12, length = 49))
#path = exp(seq(-6, -12, length = 7))
set.seed(1349)
fold = sample(rep(1:10, length = sum(subset) - length(test)))
pf = rep(NA, length(events))
names(pf) = events
ridge = naive
for (e in events) {
    print(e)
    fit = cv.glmnet(x.b, y[-test] == e, family = 'binomial',
        lambda = path, alpha = 0, standardize = FALSE, foldid = fold)
    if (fit$lambda.min == min(fit$lambda)) print(log(fit$lambda.min))
    if (fit$lambda.min == max(fit$lambda)) print(log(fit$lambda.min))
    ridge[, e] =
        predict(fit, diag(ncol(x.b)), s = 'lambda.min', type = 'response')
    pf[e] = fit$lambda.min
}
ridge = ridge/rowSums(ridge)
Sys.time() - Sys.time

# fit random effects estimator for each outcome on training data
Sys.time = Sys.time()
random = naive
for (e in events) {
    print(e)
    fit = glmer((y[-test] == e) ~ (1|batter.b),
        family = binomial(link = 'probit'))
    random[, e] = predict(fit,
        data.frame(batter.b = sort(unique(batter.b))), type = 'response')
}
random = random/rowSums(random)
Sys.time() - Sys.time

# compute test error
count.test = aggregate(y[test], by = list(batter[test]), table)
pa.test = rowSums(count.test[,-1])
y.test = count.test[,-1]/rowSums(count.test[,-1])
names(pa.test) = rownames(y.test) = sort(unique(data$BAT_ID[subset][test]))
match = names(which(pa.test > 30))
error = function(pred) {
    colSums((pred[match, ] - y.test[match, ])^2*pa.test[match])/sum(pa.test)}




# FIGURE AND TABLES FOR METHODS SECTION

# Produce Figure 3 in the paper
pdf('figs/regul-as-regre.pdf', width = 9, height = 3)
stat = '1B'
lim = c(.11, .24)
par(mfrow = c(1, 3))
plot(naive[, stat], regressed[, stat], col = 'dodgerblue',
    xlab = 'Naive estimator for 1B rate', xlim = lim,
    ylab = 'Regressed estimator for 1B rate', ylim = lim,
    main = '(a) Regressed vs. naive')
points(naive[pa < 20, stat], regressed[pa < 20, stat], col = 'darkorange',
    pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(y == stat), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean 1B rate', '< 20 PA', '> 19 PA'),
    pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
    col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
plot(naive[, stat], ridge[, stat], col = 'dodgerblue',
    xlab = 'Naive estimator for 1B rate', xlim = lim,
    ylab = 'Ridge estimator for 1B rate', ylim = lim,
    main = '(b) Ridge vs. naive')
points(naive[pa < 20, stat], ridge[pa < 20, stat], col = 'darkorange',
    pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(y == stat), lty = 2)
plot(regressed[, stat], ridge[, stat], col = 'dodgerblue',
    xlab = 'Regressed estimator for 1B rate', xlim = lim,
    ylab = 'Ridge estimator for 1B rate', ylim = lim,
    main = '(c) Ridge vs. regressed')
points(regressed[pa < 20, stat], ridge[pa < 20, stat], col = 'darkorange',
    pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(y == stat), lty = 2)
dev.off()

# Produce Table 2 in the paper
table2 = 100*sqrt(rbind(error(naive), error(regressed), error(ridge)))[, 
    c('G', 'F', 'K', 'BB', 'HBP', '1B', '2B', '3B', 'HR')]
rownames(table2) = c('Naive', 'Regressed', 'Ridge')
print(xtable(table2, digits = 2, display = c('s', rep('f', 9))),
    only.contents = TRUE, file = 'tabs/regul-as-regre.tex',
    include.colnames = FALSE, hline.after = c())

# Produce Table 3 in the paper
table3 = 100*sqrt(rbind(error(regressed), error(random)))[, 
    c('G', 'F', 'K', 'BB', 'HBP', '1B', '2B', '3B', 'HR')]
rownames(table3) = c('Regressed', 'Random')
print(xtable(table3, digits = 2, display = c('s', rep('f', 9))),
    only.contents = TRUE, file = 'tabs/regul-vs-rando.tex',
    include.colnames = FALSE, hline.after = c())




# CODE TO PRODUCE FIGURES AND TABLES IN RESULTS SECTION

weights = c(.881, 1.256, 1.594, .687, 0, 0, .718, 2.065, 0)
names(weights) = levels(y)


# VALIDATION

# Randomly select test set for holding out.
# If batter, pitcher have same handedness, put PA in test set with prob. 10%
# If batter, pitcher have oppo handedness, put PA in test set with prob. 90%
set.seed(4964)
test = which(rbinom(sum(subset), 1, .1 + .8*hand) == 1)
fold = sample(rep(1:10, length = sum(subset) - length(test)))

# Fit naive estimator of wOBA on training data
count = aggregate(y[-test], by = list(batter[-test]), table)
pa = rowSums(count[,-1])
naive = count[,-1]/rowSums(count[,-1])
names(pa) = rownames(naive) = sort(unique(data$BAT_ID[subset][-test]))
woba.naive = colSums(weights*t(naive))

# Marcel (at Tango's request)
pa.to.add = c(296, 1114, 576, 98, 0, 0, 257, 132, 59)
league.ave = c(table(y)/length(y))
count.marcel = t(t(count[,-1]) + pa.to.add*league.ave)
pa.marcel = t(t(matrix(pa,nrow(count.marcel),ncol(count.marcel))) + pa.to.add)
marcel = count.marcel/pa.marcel
rownames(marcel) = rownames(naive)
woba.marcel = colSums(weights*t(marcel))

# Fit naive estimator of wOBA for pitchers on training data (not used)
count.p = aggregate(y[-test], by = list(pitcher[-test]), table)
pa.p = rowSums(count.p[,-1])
naive.p = count.p[,-1]/rowSums(count.p[,-1])
names(pa.p) = rownames(naive.p) = sort(unique(data$PIT_ID[subset][-test]))
woba.naive.p = colSums(weights*t(naive.p))

# Fit null estimator of wOBA on training data
null = naive
for (e in events) {
    null[, e] = mean(y[-test] == e)
}

# Fit regressed estimator of wOBA on training data
regressed = naive
for (e in events) {
    priorMean = mean(y[-test] == e)
    var = priorMean*(1-priorMean)/pa
    priorVar = estimatePopulationVariance(naive[, e], var, priorMean)
    regressed[, e] = (naive[, e]/var + priorMean/priorVar)/(1/var + 1/priorVar)
}
regressed = regressed/rowSums(regressed)
woba.regressed = colSums(weights*t(regressed))

# Fit True wOBA on training data
prob.true = matrix(NA, length(test), length(events))
colnames(prob.true) = events
for (e in events) {
    fit = cv.glmnet(x[-test, ], y[-test] == e, family = 'binomial',
        lambda = path, alpha = 0, standardize = FALSE, foldid = fold)
    prob.true[, e] =
        predict(fit, x[test, ], type = 'response', s = 'lambda.min')[, 1]
}
pred.true = aggregate(prob.true, by = list(batter[test]), mean)[, -1]
rownames(pred.true) = sort(unique(data$BAT_ID[subset][test]))
woba.true = colSums(weights*t(pred.true))

# Fit linear mixed effects model for wOBA on training data
hand.t = hand[-test]
home.t = home[-test]
stadium.t = as.factor(data$HOME_TEAM_ID[subset])[-test]
batter.t = as.factor(data$BAT_ID[subset])[-test]
pitcher.t = as.factor(data$PIT_ID[subset])[-test]
fit.mixed = lmer(weights[y][-test] ~ hand.t + home.t + stadium.t +
    (1|batter.t) + (1|pitcher.t))
x.test = data.frame(hand.t = hand[test], home.t = home[test],
    stadium.t = as.factor(data$HOME_TEAM_ID[subset])[test],
    batter.t = as.factor(data$BAT_ID[subset])[test],
    pitcher.t = as.factor(data$PIT_ID[subset])[test])
match = intersect(names(pa), names(pa.test))
x.test$batter.t[!is.element(x.test$batter.t,match)]=names(which.min(pa[match]))
x.test$pitcher.t[!is.element(x.test$pitcher.t, match)] =
    names(which.min(pa.p[match]))
match = names(which(pa.test > 100))
pred.mixed = predict(fit.mixed, newdata = x.test)
woba.mixed = aggregate(pred.mixed, by = list(batter[test]), mean)[, -1]
names(woba.mixed) = names(woba.true)

# Compute test error of all methods
count.test = aggregate(y[test], by = list(batter[test]), table)
pa.test = rowSums(count.test[,-1])
y.test = count.test[,-1]/rowSums(count.test[,-1])
names(pa.test) = rownames(y.test) = sort(unique(data$BAT_ID[subset][test]))
error = function(pred) {
    colSums((pred[match, ] - y.test[match, ])^2*pa.test[match])/sum(pa.test)}
woba.test = colSums(weights*t(y.test))

error = function(pred) {
    predmatch = pred[match]
    predmatch[is.na(predmatch)] = mean(woba.regressed)
    c(wtd.mean((predmatch - woba.test[match])^2, weights = pa.test[match]),
     sqrt(wtd.var((predmatch - woba.test[match])^2, weights = pa.test[match])/length(match)))
}

error(woba.naive)
error(woba.regressed)
error(woba.marcel)
error(woba.true)
error(woba.mixed)




# True wOBA (results of fitting to entire dataset)

# Fit naive estimator of wOBA to all 2015 data
count = aggregate(y, by = list(batter), table)
pa = rowSums(count[,-1])
naive = count[,-1]/rowSums(count[,-1])
names(pa) = rownames(naive) = sort(unique(data$BAT_ID[subset]))
woba.naive = colSums(weights*t(naive))

# Fit naive estimator of wOBA against to all 2015 data
count.p = aggregate(y, by = list(pitcher), table)
pa.p = rowSums(count.p[,-1])
naive.p = count.p[,-1]/rowSums(count.p[,-1])
names(pa.p) = rownames(naive.p) = sort(unique(data$PIT_ID[subset]))
woba.naive.p = colSums(weights*t(naive.p))


# Fit True wOBA to all 2015 data
set.seed(1987)
fold = sample(rep(1:10, length = nrow(x)))
pred.true = matrix(NA, ncol(x), length(events))
colnames(pred.true) = events
for (e in events) {
    fit = cv.glmnet(x, y == e, family = 'binomial', lambda = path,
        alpha = 0, standardize = FALSE, foldid = fold)
    pred.true[, e] = predict(fit, rbind(0, 0, cbind(.5, .5, diag(ncol(x)-2))),
        type = 'response', s = 'lambda.min')[, 1]
}
rownames(pred.true) = c('hand', 'home',
    as.character(sort(unique(data$HOME_TEAM_ID[subset]))),
    as.character(sort(unique(data$BAT_ID[subset]))),
    as.character(sort(unique(data$PIT_ID[subset]))))
woba.true = colSums(weights*t(pred.true))
woba.true.bat = woba.true[sort(unique(batter))]
woba.true.pit = woba.true[sort(unique(pitcher))]

# Produce Figure 4 in the paper
pdf('figs/true-woba.pdf', height = 5, width = 10)
par(mfrow = c(1, 2))
bat.match = intersect(names(woba.naive), names(woba.true.bat))
few.pa = intersect(bat.match, names(which(pa < 20)))
lim = c(.200, .452)
axis = seq(.2, .45, length = 6)
axis.lab = c('.200', '.250', '.300', '.350', '.400', '.450')
plot(woba.naive[bat.match], woba.true.bat[bat.match], col = 'dodgerblue',
    xlim = lim, xlab = 'Observed wOBA',
    ylim = lim, ylab = 'True wOBA',
    main = '(a) True vs. Observed wOBA for batters', axes = FALSE)
axis(1, at = axis, lab = axis.lab)
axis(2, at = axis, lab = axis.lab)
points(woba.naive[few.pa], woba.true.bat[few.pa], col = 'darkorange', pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(woba.true.bat), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean True wOBA', '< 20 PA', '> 19 PA'),
    pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
    col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
pit.match = intersect(names(woba.naive.p), names(woba.true.pit))
few.pa = intersect(pit.match, names(which(pa.p < 20)))
plot(woba.naive.p[pit.match], woba.true.pit[pit.match], col = 'dodgerblue',
    xlim = lim, xlab = 'Observed wOBA against',
    ylim = lim, ylab = 'True wOBA against',
    main = '(b) True vs. Observed wOBA against pitchers', axes = FALSE)
axis(1, at = axis, lab = axis.lab)
axis(2, at = axis, lab = axis.lab)
points(woba.naive.p[few.pa], woba.true.pit[few.pa], col = 'darkorange',
    pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(woba.true.pit), lty = 2)
dev.off()

# Produce Figure 4(a) from paper for slides
pdf('figs/true-woba-batters.pdf', height = 5, width = 5)
bat.match = intersect(names(woba.naive), names(woba.true.bat))
few.pa = intersect(bat.match, names(which(pa < 20)))
lim = c(.200, .452)
axis = seq(.2, .45, length = 6)
axis.lab = c('.200', '.250', '.300', '.350', '.400', '.450')
plot(woba.naive[bat.match], woba.true.bat[bat.match], col = 'dodgerblue',
    xlim = lim, xlab = 'Observed wOBA',
    ylim = lim, ylab = 'True wOBA', axes = FALSE)
axis(1, at = axis, lab = axis.lab)
axis(2, at = axis, lab = axis.lab)
points(woba.naive[few.pa], woba.true.bat[few.pa], col = 'darkorange', pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(woba.true.bat), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean True wOBA', '< 20 PA', '> 19 PA'),
    pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
    col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
dev.off()

# Produce Figure 4(b) from paper for slides
pdf('figs/true-woba-pitchers.pdf', height = 5, width = 5)
pit.match = intersect(names(woba.naive.p), names(woba.true.pit))
few.pa = intersect(pit.match, names(which(pa.p < 20)))
plot(woba.naive.p[pit.match], woba.true.pit[pit.match], col = 'dodgerblue',
    xlim = lim, xlab = 'Observed wOBA against',
    ylim = lim, ylab = 'True wOBA against', axes = FALSE)
axis(1, at = axis, lab = axis.lab)
axis(2, at = axis, lab = axis.lab)
points(woba.naive.p[few.pa], woba.true.pit[few.pa], col = 'darkorange',
    pch = 19)
abline(0, 1, col = 'forestgreen')
abline(h = mean(woba.true.pit), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean True wOBA', '< 20 PA', '> 19 PA'),
    pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
    col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
dev.off()

# Produce Table 5 in paper
head(-sort(-woba.true.bat))
tail(-sort(-woba.true.bat))
head(sort(woba.true.pit))
tail(sort(woba.true.pit))

# Produce Table 6 in paper
qual.bat = intersect(bat.match, names(which(pa > 500)))
head(-sort(woba.naive[qual.bat] - woba.true.bat[qual.bat]))
tail(-sort(woba.naive[qual.bat] - woba.true.bat[qual.bat]))
qual.pit = intersect(pit.match, names(which(pa.p > 500)))
head(sort(woba.true.pit[qual.pit] - woba.naive.p[qual.pit]))
tail(sort(woba.true.pit[qual.pit] - woba.naive.p[qual.pit]))

