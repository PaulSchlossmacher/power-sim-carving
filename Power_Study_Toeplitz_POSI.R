#Trying to create a more extreme Toeplitz example with less noise to 
#encourage Lasso to use more variables:
# Clear all variables
rm(list = ls())


#Local, user specific path, that should work for both of us:

#save.image(file='myEnvironment_nsim200_6fraqs.RData')
#load('myEnvironment.RData')



# --------------- Create plots --------------

data_Power_POSI <- data.frame(
  Fraq=fraq.vec,
  "Avg Power Drysdale" = full_power_avg_D
)

data_Power_long_POSI <- tidyr::gather(data_Power_POSI, "Type", "Value", -Fraq)

PowerPlot_POSI<-ggplot(data_Power_long_POSI, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() +
  labs(title = "Average Power - with Bonferroni correction",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels=c('POSI'))

ggsave("PowerPlot_POSI.png", plot = PowerPlot_POSI, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")


data_TypeI_POSI <- data.frame(
  Fraq=fraq.vec,
  "Avg Type I Error rate Drysdale" = full_type1_error_avg_D
)

data_TypeI_long_POSI <- tidyr::gather(data_TypeI_POSI, "Type", "Value", -Fraq)

TypeIPlot_POSI<-ggplot(data_TypeI_long_POSI, aes(x = Fraq, y = Value, color = Type)) +
  geom_line() +
  labs(title = "Average Type I Error Rate - with Bonferroni correction",
       x = "Fraq", y = "Value") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels=c('POSI'))

ggsave("TypeIPlot_POSI.png", plot = TypeIPlot_POSI, width = 8, height = 6,
       units = "in", dpi = 300, bg = "#F0F0F0")
