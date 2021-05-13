test_that("plot_connectivity is working", {
  df <- data.frame(Median = c(10, 12),
                   Min = c(5, 6),
                   Lower = c(7, 8),
                   Upper = c(12, 13),
                   Max = c(14, 15),
                   Group = c("r_name", "t_name"))
  palette <- c("#5BBCD6", "#FF0000", "#00A08A", "#F2AD00")
  p <- df %>%
    ggplot(aes(x = Group, fill = Group)) +
    geom_boxplot(aes(ymin=Min, lower=Lower, middle=Median, upper=Upper, ymax=Max, group = Group), stat = "identity") +
    theme_minimal() +
    theme(legend.position="none") +
    ylab("Measured Connectivity") +
    ggtitle("mod_name") +
    theme(axis.text.x=element_text(angle=30))+
    scale_fill_manual(values=palette)
  expect_is(p, "gg")
})

test_that("plot_permutations is working", {
  mdc_type == "difference"
  mod_name = c("X1", "X2", "X3")
  mdc_value = 10
  df <- data.frame(X1 = sample(10,10),
                   X2 = sample(10,10),
                   X3 = sample(10,10))
  color <- ("#5BBCD6")
  p <- df %>%
    ggplot(aes(x=df)) +
    theme_minimal() +
    geom_density(color="white", fill=color) +
    ggtitle(mod_name) +
    ylab("Probability Density") +
    xlab("Permutated Differential Connectivity") +
    geom_vline(xintercept=mdc_value, linetype="dotted", size=1) +
    theme(axis.title.x=element_text(vjust=-1))
  expect_is(p, "gg")
})

