##### Script for plotting matrix representation of the feature information

library(ProductFormFA)
library(ggplot2)
library(scales)
library(ggthemes)
library(dplyr)
library(latex2exp)

set.seed(1234)

# Generate from buffet procedure from beginning
# (returns the features, the number of new features for new customers, 
# counts of observed features)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 2, theta = 5, n = 10, lambda = 30)
plot_trajectory(buff_poiss_bb$features)

# Matrix of order-of-appearance features from the buffet
ooa_mat_poiss_bb <- create_features_matrix(buff_poiss_bb$features)
mat_plot_ <- plot_binary_matrix(ooa_mat_poiss_bb)
n <- nrow(ooa_mat_poiss_bb)
K <- ncol(ooa_mat_poiss_bb)


mat_plot <- mat_plot_ + 
  scale_fill_manual(values=c("white", "black")) +
  coord_equal() + 
  labs(x= expression(italic("\u2113")), y= expression(italic("i"))) + 
  #ylab(TeX(r"($ i $)"))+
  #xlab(TeX(r"($ \ell $)")) +
  theme_light() + 
  theme(panel.grid=element_blank()) + 
  theme(panel.border=element_blank()) +
  #theme(axis.ticks =element_blank()) + 
  #theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_x_discrete(labels = 1:K, position = "top") +
  scale_y_discrete(labels = seq(n,1, by = -1)) +
  theme(legend.position="none") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "points"))


ggsave(filename = "Plots_paper/matrix_representation.pdf", width = 6, height = 3, dpi = 300, units = "in", device='pdf')


##########################################################
#### for Talk ####################
#########################################

library(ProductFormFA)
library(ggplot2)
library(scales)
library(ggthemes)
library(dplyr)
library(latex2exp)

set.seed(1234)

# Generate from buffet procedure from beginning
# (returns the features, the number of new features for new customers, 
# counts of observed features)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 2, theta = 5, n = 4, lambda = 12)
plot_trajectory(buff_poiss_bb$features)

# Matrix of order-of-appearance features from the buffet
ooa_mat_poiss_bb <- create_features_matrix(buff_poiss_bb$features)
mat_plot_ <- plot_binary_matrix(ooa_mat_poiss_bb)
n <- nrow(ooa_mat_poiss_bb)
K <- ncol(ooa_mat_poiss_bb)


mat_plot <- mat_plot_ + 
  scale_fill_manual(values=c("white", "black")) +
  coord_equal() + 
  labs(x= "s", y= "i") + 
  #ylab(TeX(r"($ i $)"))+
  #xlab(TeX(r"($ \ell $)")) +
  theme_light() + 
  theme(panel.grid=element_blank()) + 
  theme(panel.border=element_blank()) +
  #theme(axis.ticks =element_blank()) + 
  #theme(axis.text.x = element_blank()) +
  scale_x_discrete(labels = 1:K, position = "top") +
  scale_y_discrete(labels = seq(n,1, by = -1)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text=element_text(size=10),
        axis.title.x.top = element_text(margin = margin(1, 0, 10, 0)),
        axis.title.y.left = element_text(margin = margin(0, 10, 0, 1))) + 
  theme(legend.position="none") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "points"))


ggsave(filename = "Plots_paper/matrix_representation_talk.pdf", width = 4, height = 2, dpi = 300, units = "in", device='pdf')

