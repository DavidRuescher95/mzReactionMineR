object <- data_normalized
assay <- "normalized"
sample_column <- "LC_MS_Nr"
design_column <- "Vector"
controls <- c("pAGT1009")
thrsh <- 0.8
id_col = "id"
rt_col = "rt"
mz_col = "mz"

detectNewCompounds <- function(
  object,
  assay,
  sample_column,
  design_column,
  controls,
  thrsh,
  id_col = "id",
  rt_col = "rt",
  mz_col = "mz"
)


tmp <- cbind(
  as.data.frame(rowData(object))[,c(id_col,rt_col,mz_col)],
  as.data.frame(assays(object)[[assay]]) %>% 
    replace(is.na(.),0) %>%
    replace(. > 0, 1)
) %>%
  pivot_longer(
    cols = -c(!!sym(id_col), !!sym(rt_col), !!sym(mz_col)),
    names_to = sample_column,
    values_to = "Value"
  ) %>%
  inner_join(
    as.data.frame(colData(object)) 
  )

# all features that are not in control

candidates_1 <- tmp %>%
  filter(
    !!sym(design_column) %in% controls & Value == 0
  ) %>%
  select(!!sym(id_col), !!sym(rt_col), !!sym(mz_col)) %>%
  unique() %>%
  mutate(
    Type = "not_in_control"
  ) %>%
  as.data.frame()

# all features not in control, but in 80% of any group

candidates_2 <- tmp %>%
  filter(
    !(!!sym(design_column) %in% controls) & Value != 0,
    !!sym(id_col) %in% candidates_1[,id_col]
  ) %>%
  group_by(!!sym(id_col), !!sym(rt_col), !!sym(mz_col), !!sym(design_column)) %>%
  summarize(
    Value = sum(Value),
    cutoff = floor(n()*thrsh)
  ) %>%
  filter(
    id %in% candidates_1$id,
    Value >= cutoff
  ) %>%
  select(!!sym(id_col), !!sym(rt_col), !!sym(mz_col), !!sym(design_column)) %>%
  unique() %>%
  mutate(
    Type = "in_vector"
  )

# filtered for relevant mz range

candidates_3 <- tmp %>%
  filter(
    id %in% candidates_2$id,
    mz > 250 & mz < 600
  ) %>%
  select(id,rt,mz) %>%
  unique() %>%
  mutate(
    Type = "mz_range"
  )

candidates <- rbind(
  candidates_1,
  candidates_2,
  candidates_3
)

rm(tmp, candidates_1, candidates_2, candidates_3)



Method3 <- rbind(
  rowData(object) %>%  %>% as.data.frame() %>% mutate(Type = "Aligned"),
  candidates
)

plot_data <- Method3 %>%
  group_by(Type) %>%
  summarize(
    n = n()
  )

ggplot(
  data = plot_data,
  aes(
    x = fct_reorder(Type,n,.desc=TRUE),
    y = n,
    fill = n
  )
) +
  geom_col(
    width = 0.5,
    color = "black"
  )  +
  geom_text(
    aes(
      label = n,
      y = n + max(n)*0.05
    ),
    size = 2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_viridis_c(option="plasma") +
  labs(
    x = "Filter stage",
    y = "Number of features"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.background = element_rect(fill = NA, color = "black"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(color = NA, fill = NA),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )
ggsave(
  paste0("Method3.png"),
  path = file.path(results_dir,"strategy/figures"),
  units = "cm",
  width = 10,
  height = 8,
  dpi = 600
)
