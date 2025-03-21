# AUTHOR: DANIEL HERNÁNDEZ CARRSCO
# DATE: 2025-03-21

# Load packages
library(ggplot2)
library(cowplot)
library(showtext)
library(latex2exp)

# Plotting preferences
theme_set(theme_bw())
theme_update(text = element_text(size = 9),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             legend.key.size = unit(0.2, "in"),
             axis.text = element_text(),
             axis.title = element_text(),
             panel.background = element_rect(fill = "white",colour = "black"), 
             strip.background = element_rect(color = "transparent", 
                                             fill = "transparent", size = 1,
                                             linetype = "solid"),
             legend.background = element_rect(fill = "transparent"),
             plot.background = element_rect(fill = "transparent", color = NA),
             strip.text = element_text( color = "black", face = "italic", size = 12),
             axis.ticks.length = unit(0.1, "cm"), 
             axis.ticks = element_line(size = 0.3, colour = "black"),
             axis.text.x = element_text(margin = unit(c(0.2,0.2,0.2,0.2), "cm"), colour = "black"), 
             axis.text.y = element_text(margin = unit(c(0.2,0.2,0.2,0.2), "cm"), colour = "black")
)

##### 1 | EFFECT OF NON-LINEAR RESPONSE ON AVERAGE LONG-TERM VITAL RATES ####

# Environmental response
vital_rate <- Vectorize(function(r,n, s) (r^s)/(n^s + r^s)) # Environmental response function (e.g., survival)
vital_rate_var <- Vectorize(function(r,s,n,amp,period, t){ # Average vital rate under seasonal fluctuations (average adter full period)
  r_t = r + amp * sin((2*pi*(1:t))/period)
  return(mean(vital_rate(r = r_t, s = s, n = n)))
})

# Second derivative of response curve
ddvital_rate <- function(r, n, s) {
  ( s * n^s * r^(s - 2) * ((s - 1) * (n^s + r^s) - 2 * s * r^s) ) / (n^s + r^s)^3
}

# Average environment (e.g., resource availability
r <- seq(from = 0, to = 5, length.out = 10000) 

# Amplitude of seasonal fluctuations
amp = 1

# Parameters determining the shape of the response curve
s = 5
n = 2

# Test  
vital_rate(r = 1:2, s = 5, n = 2)
vital_rate_var(r = 1:2, s = 5, n = 2, amp = 1, period = 100, t  = 100)
ddvital_rate(r = 1:2, s = 5, n = 2)

# Create data frame
functions.data  <- data.frame(r = r, # Average environment
                              vital_rate = vital_rate(r,s = s, n = n), # vital rate under average conditions
                              vital_rate_var = vital_rate_var(r = r,s = s,n = n, amp = amp, t = 100, period = 100), # average vital rate under seasonal fluctuations
                              ddvital_rate = ddvital_rate(r,s = s, n = n))


# Plot second derivative
plot2d <- ggplotGrob(ggplot(functions.data, aes(x = r, y = ddvital_rate)) + 
                       theme_void() +
                       geom_area(data = subset(functions.data, ddvital_rate < 0),
                                 col = "grey20", 
                                 fill = '#2a737a', 
                                 size = 0.3, 
                                 alpha = 0.5) + 
                       geom_area(data = subset(functions.data, ddvital_rate > 0),
                                 col = "grey20", 
                                 fill = 'darkorange', 
                                 size = 0.3,
                                 alpha = 0.5) + 
                       ylab(TeX("$\\lambda''(x)$")) +
                       geom_hline(yintercept = 0, size = 0.3) +
                       xlab(TeX("$x$")) +
                       scale_x_continuous(expand = c(0,0)) +
                       theme(legend.title = element_blank(),
                             aspect.ratio = 1, 
                             panel.grid = element_blank(),
                             plot.title = element_text(size = 7, hjust = 0.5),
                             axis.ticks = element_line(size = 0.2),
                             axis.title.x = element_text(size = 6, margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")),
                             axis.title.y = element_text(size = 6, margin = unit(c(0.1,0.1,0.1,0.1), units = "cm"), angle = 90),
                             axis.ticks.length = unit(0.1, "cm"), 
                             panel.border = element_rect(fill = "transparent"),
                             panel.background = element_rect(fill = "grey98"),
                             plot.margin = unit(c(0,0,0,0), units = "cm")) +
                       annotate(geom = "text", x = 2, y = 0.3, 
                                label = TeX("$\\lambda''(x)>0$", output = "character"), 
                                hjust = "left",
                                size = 2,
                                parse = TRUE) + 
                       annotate(geom = "text", x = 3.5, y = -0.3, 
                                label = TeX("$\\lambda''(x)<0$", output = "character"), 
                                hjust = "left",
                                size = 2,
                                parse = TRUE) + 
                       ggtitle("Convexity"))

# Plot average vital rate
(plot_jensens <-
  ggplot(functions.data, aes(x = r, y = vital_rate)) +
  geom_line(aes(lty = "Constant"), size = 0.3) + 
  geom_line(data = subset(functions.data, ddvital_rate > 0),
            aes(y = vital_rate_var, lty = "Seasonal"),
            size = 0.3, 
            colour = 'darkorange') +
  geom_line(data = subset(functions.data, ddvital_rate < 0),
            aes(y = vital_rate_var, lty = "Seasonal"),
            size = 0.3, colour = '#2a737a') +
  scale_linetype_manual(values = c(1,5)) +
  ylab(TeX("Average demographic rate  ( $\\bar{\\lambda}$)")) +
  xlab(TeX("Environment (x)")) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(hjust = 0),
        legend.position = c(0.2,0.8),
        legend.direction = "vertical",
        panel.border = element_blank(),
        aspect.ratio = 1) + 
  coord_cartesian(ylim = c(0, 1)) +
  guides(colour = 'none' ) +
  annotation_custom(grob = plot2d, xmin = 2, xmax = 4.8, ymin = -0.05, ymax = 0.6))



#### 2 | EFFECT OF FLUCTUATING GROWTH-RATES ON LONG TERM POPULATION GROWTH ####

set.seed(1)

# Parameters
base_lambda = 1.01 # Average lambda
amp_seas = 0.2 # Amplitude of seasonal fluctuations
amp_noise = 0.02
period = 12 # Period of fluctuations
t = period*10 # Total time
N0 = 10 # Initial abundance


# Growth rates
lambda_0 = base_lambda + amp_noise*rnorm(t)
lambda_S = lambda_0 + amp_seas*cos((2*pi*(1:t)/period))

# check average growth rates of both populations is the same
mean(lambda_0)
mean(lambda_S)

# Initialize population abundance
seas <- c(N0); N_s = N0
const <- c(N0); N_c = N0

# Apply population growth
for (i in lambda_0) {
  N_c = N_c*i
  const = c(const, N_c)
}

for (i in lambda_S) {
  N_s = N_s*i
  seas = c(seas, N_s)
}

# Create data frame
growth.data = rbind(data.frame(N = const, 
                               G = c(lambda_0,NA),
                               type = 'constant',
                               time = 0:t),
                    data.frame(N = seas, 
                               type = 'seasonal',
                               G = c(lambda_S, NA),
                               time = 0:t))

# Abundance
gA <- ggplot(growth.data, aes(x = time/12, y = N)) +
    theme(legend.position = '',
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    geom_line(aes(lty = type), size = 0.3) +
    scale_linetype_manual( values = c(1,5)) +
    scale_y_continuous('Population size') + 
    scale_x_continuous(expand = c(0,0), breaks = 0:100)

# Growth rate
gB <- ggplot(growth.data, aes(x = time/12, y = G)) +
    geom_hline(yintercept = 1, lty = 1, colour = 'grey70') +
    geom_line(aes(lty = type), size = 0.3) +
    theme(legend.position = '',
          panel.grid = element_blank()) +
    scale_y_continuous("Growth rate", breaks = c(0.8, 1,1.2)) +
    scale_linetype_manual( values = c(1,5)) +
    scale_x_continuous('Time (years)', breaks = 0:100, expand = c(0,0))

gAB <- plot_grid(gA, gB, ncol = 1, align = 'v', rel_heights = c(1,0.67))


# Combine plots and save
plot_grid(plot_jensens, gAB, labels = c("A", "B"), ncol = 2, rel_widths =  c(1, 1.3), label_size = 10)
ggsave("./demographic_rates.pdf", width = 19, height = 8, units = 'cm', dpi = 600)



