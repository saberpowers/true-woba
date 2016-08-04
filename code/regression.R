# wOBA Regression Project
# Eli Shayer
# regression.R
# Regress wOBA to the mean for 2015 batters

require(xtable)

# set the working directory
# Assumes that files have been processed with Chadwick
# as described in Analyzing Baseball Data with R (appendix 1)

##### Constants #####
# the year of the analysis
# currently only 2015, but in this form for easier extension
year = 2015
minPa = 300
paQual = 10

##### Data Reading #####
# read the play by play data for the current year
# headers come from http://chadwick.sourceforge.net/doc/cwevent.html
pbp = read.csv(paste0(c("data/all", year, ".csv"), collapse=""),
               header = F, stringsAsFactors = F)
colnames(pbp) = read.csv('data/fields.csv')$Header

# remove all PA in which a pitcher is the hitter
pbp = pbp[pbp$BAT_FLD_CD != 1,]

##### Data Processing #####
# event type vectors for Retrosheet data
# http://chadwick.sourceforge.net/doc/cwevent.html
pa = c(2, 3, 14:23)
wOBADem = c(2, 3, 14, 16:23)
bip =  c(2, 18:22)
hip = c(20:22)


# get the list of qualifying batters
batters = aggregate(pbp$EVENT_CD %in% pa, by = list(pbp$BAT_ID), sum)
batters = batters[batters$x >= paQual,]
allBattersPA = batters$x
allBatters = batters$Group.1
allQual = is.element(pbp$BAT_ID, allBatters)

# threshold batters based on minimum plate appearances, store ids
batters = batters[batters$x >= minPa,]$Group.1

# helper function to get the count of a logical vector by batter
# with the batter vector determined by the passed in value
# logicalVector is the lenght of pbp with TRUE if to be counted
countByBatter = function(logicalVector, players) {
  data = aggregate(logicalVector, by = list(pbp$BAT_ID), sum)
  return(data[data$Group.1 %in% players,]$x)
}

# A function to get the league variance of a binary stat by iteratively
# combining the estimates of league variance from individual players
# until the estimate converges.
getStatVariance = function(denX, numX) {
  # a data frame with seven columns as follows:
  # player id, denominator count, numerator count, stat rate,
  # estimate of pop variance based on player,
  # variance of the estimate of pop variance based on player,
  # the variation in the player's skill.
  # initialize to 0 for an empty sum for counting purposes
  df = as.data.frame(matrix(0, length(batters), 7))
  colnames(df) = c("player", "den", "num", "stat",
                   "varEstimate", "varVarEstimate", "varSkill")
  
  # store the player ids in the first column of the data frame
  df$player = batters
  
  # specific to strikeouts and pa, reduced to just the qualifying batters
  df$num = countByBatter(numX & denX, batters)
  df$den = countByBatter(denX, batters)
  
  # generate and store the stat rate
  df$stat = df$num / df$den
  
  # store the league average skill and number of events in sample space
  statBar = sum(df$num) / sum(df$den)
  nTotal = sum(df$den)
  
  # All of the following comes from The Book page 387
  # estimate of the league variation for each player
  df$varEstimate = (df$stat - statBar) ^ 2 -
    (1 - df$den / nTotal) * df$stat * (1 - df$stat) / df$den
  
  # initialize the variation in the skill to one, to be found iteratively
  leagueVar = 1
  
  # iteratively find the league variation, arbitrary cap on iterations
  iterationsCount = 0
  while (iterationsCount < 1e5) {
    prevLeagueVar = leagueVar
    
    # variance in the estimate of league variance for each player
    df$varVarEstimate = sqrt(2) * (df$stat * (1 - df$stat) / df$den + leagueVar)
    
    leagueVar = sum(df$varEstimate / df$varVarEstimate ^ 2) /
      sum(1 / df$varVarEstimate ^ 2)
    
    # arbitrary percent change threshold at which to exit the loop
    if (abs((prevLeagueVar - leagueVar) / leagueVar) < 1e-15) {
      iterationsCount = 1e5
    }
    
    iterationsCount = iterationsCount + 1
  }
  
  return(leagueVar)
}

# a function to get a statistic regressed to the mean
regressToMean = function(denX, numX, batters, qual = TRUE, num = 150) {
  # if the event is not possible (no one did it) return all zeroes
  if (sum(denX & numX) == 0) {
    print('no events in sample space')
    if (length(batters) == 1) {
      return(as.vector(matrix(0, length(allBatters), 1)))
    } else {
      return(as.vector(matrix(0, length(batters), 1)))
    }
  }
  
  leagueVar = getStatVariance(denX, numX)
  statBar = sum(numX) / sum(denX)
  
  if (length(batters) == 1) {
    result = as.data.frame(matrix(0, length(allBatters), 4))
    colnames(result) = c("player", "n", "stat", "varSkill")
    result$player = allBatters
    result$n = countByBatter(denX, allBatters)
    result$stat = countByBatter(numX & denX, allBatters) / result$n
  } else {
    result = as.data.frame(matrix(0, length(batters), 4))
    colnames(result) = c("player", "n", "stat", "varSkill")
    result$player = batters
    result$n = 150

    stat = aggregate(numX[qual&denX], by=list(pbp$BAT_ID[qual&denX]),
                     function(x) {sum(x[1:num])})$x/num
    result$stat = ifelse(is.na(stat), 0, stat)
  }
  
  # get the variance in the player's skill
  result$varSkill = statBar * (1 - statBar) / result$n
  
  # return the stat regressed to the mean, with NaN replaced as 0, as a vector
  return(unlist(rapply(list(((statBar / leagueVar) + (result$stat / result$varSkill)) /
    ((1 / leagueVar) + (1 / result$varSkill))), f=function(x) {
      ifelse(is.nan(x), 0, x)
    }, how="replace")))
}

# helper functions for producing logical vectors for the pbp data
isEventType = function(type) {
  return(pbp$EVENT_CD %in% type)
}
isBattedBallType = function(type) {
  return(pbp$BATTEDBALL_CD %in% type)
}

# batted ball types
battedBallTypes = c('F', 'L', 'G', 'P')

datawOBA = as.data.frame(matrix(0, 6, 3))
colnames(datawOBA) = c('weight', 'isInPlay', 'num')
rownames(datawOBA) = c('ubb', 'hbp', 'single', 'double', 'triple', 'hr')

# 2015 wOBA weights: ubb, hbp,  1b,   2b,     3b,     hr
datawOBA$weight = c(.687,	.718,	.881,	1.256,	1.594,	2.065)

# set whether each wOBA component involves a ball in play
datawOBA$isInPlay = c(rep(F, 2), rep(T, 4))

# set the numerator to the appropriate event codes for all wOBA components
datawOBA$num = c(14, 16, 20, 21, 22, 23)

# calculate regressed to the mean wOBA and actual wOBA
wOBA = 0
wOBAReg = 0
# simple regression to the mean for each component
for (i in 1:nrow(datawOBA)) {
  event = is.element(pbp$EVENT_CD, datawOBA$num[i])
  wOBA = wOBA + datawOBA$weight[i] * countByBatter(isEventType(datawOBA$num[i]), allBatters)
  wOBAReg = wOBAReg + datawOBA$weight[i] * regressToMean(paX, event, FALSE, TRUE)
}

# divide wOBA by the denominator for each batter
wOBA = wOBA / countByBatter(isEventType(wOBADem), allBatters)

# get BABIP
BABIP = countByBatter(hip, allBatters) /
  countByBatter(bip, allBatters)

# wOBA regression to the mean in a data frame
result = as.data.frame(allBatters)
result = cbind(result, allBattersPA, wOBA, wOBAReg, BABIP)
colnames(result) = c('player', 'PA', 'wOBA', 'estimatedwOBA', 'BABIP')

##### TABLE 1
paX = is.element(pbp$EVENT_CD, c(2, 3, 14:23))
onbase = is.element(pbp$EVENT_CD, c(14:17, 20:23))
bip = is.element(pbp$EVENT_CD, c(2, 18:22))
hip = is.element(pbp$EVENT_CD, 20:22)
BIP = aggregate(bip, by = list(pbp$BAT_ID), sum)
minBatters = BIP$Group.1[BIP$x > 300]
qual = is.element(pbp$BAT_ID, minBatters)

first150 = function(x) {sum(x[1:150])}
second150 = function(x) {sum(x[151:300])}

table1 = as.data.frame(matrix(NA, 7, 3))
colnames(table1) = c('Metric', 'RMSE', 'Correlation')

# babip
babip1 = aggregate(hip[qual&bip], by = list(pbp$BAT_ID[qual&bip]), first150)$x/150
babip2 = aggregate(hip[qual&bip], by = list(pbp$BAT_ID[qual&bip]), second150)$x/150
babipRttm = regressToMean(bip, hip, minBatters, qual)

table1[1,] = c('BABIP Naive', sqrt(mean((babip2 - babip1) ^ 2)), cor(babip2, babip1))
table1[2,] = c('BABIP RTTM', sqrt(mean((babip2 - babipRttm) ^ 2)), cor(babip2, babipRttm))

#obp
obp1 = aggregate(onbase[qual&paX], by = list(pbp$BAT_ID[qual&paX]), first150)$x/150
obp2 = aggregate(onbase[qual&paX], by = list(pbp$BAT_ID[qual&paX]), second150)$x/150
obpRttm = regressToMean(paX, onbase, minBatters, qual)

table1[3,] = c('OBP Naive', sqrt(mean((obp2 - obp1) ^ 2)), cor(obp2, obp1))
table1[4,] = c('OBP RTTM', sqrt(mean((obp2 - obpRttm) ^ 2)), cor(obp2, obpRttm))

#woba
wOBA1 = 0
wOBA2 = 0
wOBAEvent = 0
wOBAEventAndType = 0

# calculate wOBA
for (i in 1:nrow(datawOBA)) {
  event = is.element(pbp$EVENT_CD, datawOBA$num[i])
  wOBA1 = wOBA1 + datawOBA$weight[i] * aggregate(event[qual&paX], by = list(pbp$BAT_ID[qual&paX]), first150)$x/150
  wOBA2 = wOBA2 + datawOBA$weight[i] * aggregate(event[qual&paX], by = list(pbp$BAT_ID[qual&paX]), second150)$x/150
  wOBAEvent = wOBAEvent + datawOBA$weight[i] * regressToMean(paX, event, minBatters, qual)
  
  if (datawOBA$isInPlay[i]) {
    for (type in battedBallTypes) {
      battedBallEvent = is.element(pbp$BATTEDBALL_CD, type)
      wOBAEventAndType = wOBAEventAndType + datawOBA$weight[i] *
        regressToMean(paX, battedBallEvent & event, minBatters, qual)
    }
  } else {
    wOBAEventAndType = wOBAEventAndType + datawOBA$weight[i] * regressToMean(paX, event, minBatters, qual)
  }
}

table1[5,] = c('wOBA Naive', sqrt(mean((wOBA2 - wOBA1) ^ 2)), cor(wOBA2, wOBA1))
table1[6,] = c('wOBA RTTM Events', sqrt(mean((wOBA2 - wOBAEvent) ^ 2)), cor(wOBA2, wOBAEvent))
table1[7,] = c('wOBA RTTM Events and Type', sqrt(mean((wOBA2 - wOBAEventAndType) ^ 2)), cor(wOBA2, wOBAEventAndType))

print(xtable(table1, digits = 3, display = c('s', rep('f', 3))),
      only.contents = TRUE, file = 'C:/Users/Eli/Dropbox/volatility/tabs/naive-vs-rttm.tex',
      include.colnames = FALSE, hline.after = c())

##### TABLE 2
PA = aggregate(paX, by = list(pbp$BAT_ID), sum)
tableTwoBatters = PA$Group.1[PA$x > 400]
qualified = is.element(pbp$BAT_ID, tableTwoBatters)

first200 = function(x) {sum(x[1:200])}
second200 = function(x) {sum(x[201:400])}

table2 = as.data.frame(matrix(NA, 3, 9))
colnames(table2) = c('G', 'F', 'K', 'BB', 'HBP', '1B', '2B', '3B', 'HR')
rownames(table2) = c('Naive error', 'Regressed error', 'Population variance')

table2events = as.data.frame(matrix(NA, nrow(pbp), 9))
colnames(table2events) = c('G', 'F', 'K', 'BB', 'HBP', '1B', '2B', '3B', 'HR')

y = as.character(pbp$EVENT_CD)

table2events$G = (y == 2 | y == 18 | y == 19) & pbp$BATTEDBALL_CD != 'G'
table2events$F = (y == 2 | y == 18 | y == 19) & pbp$BATTEDBALL_CD == 'G' & !pbp$SH_FL
table2events$K = y == 3
table2events$BB = y == 14
table2events$HBP = y == 16
table2events[,'1B'] = y == 20
table2events[,'2B'] = y == 21
table2events[,'3B'] = y == 22
table2events[,'HR'] = y == 23

for (i in 1:ncol(table2events)) {
  stat1 = aggregate(table2events[,i][qualified&paX], by = list(pbp$BAT_ID[qualified&paX]), first200)$x/200
  stat2 = aggregate(table2events[,i][qualified&paX], by = list(pbp$BAT_ID[qualified&paX]), second200)$x/200
  statRttm = regressToMean(paX, table2events[,i], tableTwoBatters, qual = qualified, num = 200)
  popVar = getStatVariance(paX, table2events[,i])
  
  table2[,i] = 100 * c(sqrt(mean((stat2 - stat1) ^ 2)), sqrt(mean((stat2 - statRttm) ^ 2)), popVar * 100)
}

print(xtable(table2, digits = 3, display = c('s', rep('f', 9))),
      only.contents = TRUE, file = 'C:/Users/Eli/Dropbox/volatility/tabs/rttm-by-event.tex',
      include.colnames = TRUE, hline.after = c())

# population variances
table3 = as.data.frame(matrix(NA, 2, 2))
colnames(table3) = c('Metric', 'Population Variance')

table3[1,] = c('BABIP', getStatVariance(bip, hip))
table3[2,] = c('OBP', getStatVariance(paX, onbase))

print(xtable(table3, digits = 3, display = c('s', rep('f', 2))),
      only.contents = TRUE, file = 'C:/Users/Eli/Dropbox/volatility/tabs/pop-var-babip-obp.tex',
      include.colnames = TRUE, hline.after = c())

##### PLOTING
# get the boundary between the first and second quartiles
firstQuartile = quantile(result$PA)[2]

##### PLOTTING
# get the boundary between the first and second quartiles
firstQuartile = quantile(result$PA)[2]

# plot observed v estimated wOBA with a 1:1 red line
pdf('figs/woba-observed-v-estimated.pdf', height = 7*(5/6), width = 7)
with(result, plot(wOBA, estimatedwOBA, axes = FALSE,
                  col = ifelse(PA >= 50, "dodgerblue", 'darkorange'),
                  pch = ifelse(PA >= 50, 1, 19),
                  xlab = 'Observed wOBA',
                  ylab = 'Estimated wOBA',
                  xlim = c(.150, .400),
                  ylim = c(.200, .400)))
axis(1, at = c(.150, .200, .250, .300, .350, .400),
     label = c('.150', '.200', '.250', '.300', '.350', '.400'))
axis(2, at = c(.200, .250, .300, .350, .400),
     label = c('.200', '.250', '.300', '.350', '.400'))
abline(0, 1, col = 'forestgreen')
abline(h = mean(result$estimatedwOBA), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean Estimated wOBA', '< 50 PA', '> 49 PA'),
       pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
       col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
dev.off()

pdf('figs/woba-observed-v-estimated-slides.pdf', height = 5, width = 6)
with(result, plot(wOBA, estimatedwOBA, axes = FALSE,
                  col = ifelse(PA >= 50, "dodgerblue", 'darkorange'),
                  pch = ifelse(PA >= 50, 1, 19),
                  xlab = 'Observed wOBA',
                  ylab = 'Estimated wOBA',
                  xlim = c(.150, .400),
                  ylim = c(.200, .400)))
axis(1, at = c(.150, .200, .250, .300, .350, .400),
     label = c('.150', '.200', '.250', '.300', '.350', '.400'))
axis(2, at = c(.200, .250, .300, .350, .400),
     label = c('.200', '.250', '.300', '.350', '.400'))
abline(0, 1, col = 'forestgreen')
abline(h = mean(result$estimatedwOBA), lty = 2)
legend('topleft', c('Diagonal y = x', 'Mean Estimated wOBA', '< 50 PA', '> 49 PA'),
       pch = c(NA, NA, 19, 1), lty = c(1, 2, NA, NA), box.lwd = -1,
       col = c('forestgreen', 'black', 'darkorange', 'dodgerblue'))
dev.off()


# plot (estimated wOBA - observed wOBA) v observed wOBA
# with a horizontal line at error = 0

# For the line of code below to work, you need to have a folder labeled figs
# in your working directory. If you uncomment the line below you also want to
# uncomment the dev.off() line below the plotting code.
pdf('figs/woba-v-babip.pdf')
with(result, plot(BABIP, estimatedwOBA - wOBA, axes = FALSE,
                  col = "dodgerblue",
                  xlab = 'Observed BABIP', xlim = c(.17, .43),
                  ylab = 'Projected increase in wOBA', ylim = c(-.11, .15)))
axis(1, at = c(.2, .25, .3, .35, .4),
     label = c('.200', '.250', '.300', '.350', '.400'))
axis(2, at = c(-.1, -.05, 0, .05, .1),
     label = c('-.100', '-.050', '0', '+.050', '+.100'))
abline(0, 0)
dev.off()

pdf('figs/woba-v-babip-slides.pdf', height = 5, width = 5)
with(result, plot(BABIP, estimatedwOBA - wOBA, axes = FALSE,
                  col = "dodgerblue",
                  xlab = 'Observed BABIP', xlim = c(.17, .43),
                  ylab = 'Projected increase in wOBA', ylim = c(-.11, .15)))
axis(1, at = c(.2, .25, .3, .35, .4),
     label = c('.200', '.250', '.300', '.350', '.400'))
axis(2, at = c(-.1, -.05, 0, .05, .1),
     label = c('-.100', '-.050', '0', '+.050', '+.100'))
abline(0, 0)
dev.off()
