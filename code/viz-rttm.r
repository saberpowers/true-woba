# code to produce pretty normal curves on slides 5-7
# and function plots on slides 14, 17

require(ggplot2)

blank =
    theme(axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    legend.position="none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank())

x = seq(-3, 6, length = 1000)
y = c(dnorm(x), dnorm(x, 2, sqrt(2)), dnorm(x, 2/3, sqrt(2/3)))
dt = data.table(rep(x, 3), y)

pdf('figs/talent.pdf', width = 28)
ggplot(dt) +
    geom_ribbon(data = dt[1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'dodgerblue', alpha = 0.5) + blank +
    coord_cartesian(ylim = range(dt$y))
dev.off()

pdf('figs/talent-skill.pdf', width = 28)
ggplot(dt) +
    geom_ribbon(data = dt[1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'dodgerblue', alpha = 0.5) +
    geom_ribbon(data = dt[1000 + 1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'darkorange', alpha = 0.5) + blank +
    coord_cartesian(ylim = range(dt$y))
dev.off()

pdf('figs/talent-skill-post.pdf', width = 28)
ggplot(dt) +
    geom_ribbon(data = dt[1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'dodgerblue', alpha = 0.5) +
    geom_ribbon(data = dt[1000 + 1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'darkorange', alpha = 0.5) +
    geom_ribbon(data = dt[2000 + 1:1000, ], aes(x = x, ymax = y), ymin = 0,
    fill = 'forestgreen', alpha = 0.5) + blank +
    coord_cartesian(ylim = range(dt$y))
dev.off()

x = seq(-6, 6, length = 100)

pdf('figs/logistic.pdf', height = 3)
plot(x, exp(x)/(1+exp(x)), type = 'l', col = 'forestgreen', lwd = 6,
    xlab = '', ylab = '', axes = FALSE)
axis(1)
axis(2, at = c(0, 0.5, 1), pos = 0)
dev.off()

pdf('figs/normalcdf.pdf', height = 3)
plot(x, pnorm(x), type = 'l', col = 'forestgreen', lwd = 6,
    xlab = '', ylab = '', axes = FALSE)
axis(1)
axis(2, at = c(0, 0.5, 1), pos = 0)
dev.off()
