# run after data processing in mainv2.R

# longitudinal trend plot, test
source('report_functions.R')
library(ggpubr)
cur_sample <- 'ISG-34'
healthy_range <- (ISG %>% filter(Group=='Healthy control'))$geomean %>% quantile(c(0.01,0.99))
df <- ISG %>% 
  prepare_timeline_table(current_sample = cur_sample, 
                         sample_info = sample_info_all %>% filter(`Patient ID` == cur_sample))
start_day <- df$`Date of sampling`[1]
df$timepoints <- difftime(df$`Date of sampling`, start_day, units="days") %>% as.double()
# ISG-1
# vertical_line_label <- c('2023-03-09', '2023-04-06', '2023-05-08', '2023-06-05', '2023-07-04', '2023-07-27', '2023-08-16', '2023-09-11', '2023-10-10', '2023-11-07', '2023-12-07', '2024-01-08')
# ISG-6
# vertical_line_label <- c('2023-03-02', '2023-03-27', '2023-04-25', '2023-05-23', '2023-06-16', '2023-06-29', '2023-07-25', '2023-08-23', '2023-09-18', '2023-10-16', '2023-11-13', '2023-12-11', '2024-01-08')
# ISG-34
vertical_line_label <- c('2023-05-10', '2023-06-08', '2023-07-06', '2023-08-05', '2023-09-06', '2023-10-31')

vertical_line <- difftime(as.Date(vertical_line_label), start_day, units = "days")

set.seed(42)
x_min = min(df$timepoints, na.rm = T) -100
x_max = max(df$timepoints, na.rm = T) +100
ggplot(df, aes(x=timepoints, y=geomean)) + 
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
            fill = "lightblue", alpha=0.05) +
  geom_point() + 
  geom_path(data = . %>% filter(connection_group), aes(group=connection_group), linetype = 2) +
  geom_vline(xintercept = vertical_line) + 
  annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
  scale_x_continuous(breaks = unique(df$timepoints), name = 'Days') + 
  geom_text_repel(aes(label=label), size = 2.5) +
  # scale_y_continuous(trans='log1p') + 
  annotate(geom="text", x=x_min, y=healthy_range[2], hjust = 0, vjust = 1, label="Range of healthy controls",
           color="blue")+
  theme_bw() +
  theme(aspect.ratio=3/15, 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank())
ggsave('figures/ISG-34 trendline.pdf')

# visualize

df <- left_join(ISG %>% select(Sample, geomean),
                NFkb %>% select(Sample, geomean), by = 'Sample') %>% 
  mutate(geomean = geomean.y/geomean.x)

NFkb_ratio <- NFkb %>% select(-geomean) %>% left_join(df, by='Sample') %>%
  filter(!(Group %in% c('Patients', 'Patient IFN stimulation')) | 
           `Subject ID` == cur_sample) %>%
  mutate(label = if_else(`Subject ID` == cur_sample, 
                         paste(`Subject ID`, Visit, sep='_'), 
                         NA))

# score_scatter_plot(NFkb_ratio) + 
#   scale_y_continuous() +
#   ylab('NFkB to IFN score ratio')
# ggsave('figures/ISG-6_NFkB_ratio scatter.pdf')

healthy_range <- (NFkb_ratio %>% filter(Group=='Healthy control'))$geomean %>% quantile(c(0.01,0.99))

df <- NFkb_ratio %>% 
  prepare_timeline_table(current_sample = cur_sample, 
                         sample_info = sample_info_all %>% filter(`Patient ID` == cur_sample))

df$timepoints <- difftime(df$`Date of sampling`, start_day, units="days") %>% as.double()
vertical_line <- difftime(as.Date(vertical_line_label), start_day, units = "days")

set.seed(42)
x_min = min(df$timepoints, na.rm = T) -100
x_max = max(df$timepoints, na.rm = T) +100
ggplot(df, aes(x=timepoints, y=geomean)) + 
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
            fill = "lightblue", alpha=0.05) +
  geom_point() + 
  geom_path(data = . %>% filter(connection_group), aes(group=connection_group), linetype = 2) +
  geom_vline(xintercept = vertical_line) + 
  annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
  scale_x_continuous(breaks = unique(df$timepoints), name = 'Days') + 
  geom_text_repel(aes(label=label), size = 2.5) +
  # scale_y_continuous(trans='log1p') + 
  annotate(geom="text", x=x_min, y=healthy_range[2], hjust = 0, vjust = 1, label="Range of healthy controls",
           color="blue")+
  theme_bw() +
  theme(aspect.ratio=3/15, 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank()) +
  ylab('NFkB to IFN score ratio')

ggsave('figures/ISG-34 NFkB_ratio trendline.pdf')

# differentially abundant ISGs
# cur_sample <- 'ISG-6'
lapply(c('ISG-6', 'ISG-34'), function(cur_sample){
  columns <- na.omit(unlist(panel_info[, 'NF-kB score'], use.names = F))
  df <- dat %>% 
    select(all_of(columns), `Subject ID`, Visit) %>%
    filter(`Subject ID` == cur_sample) %>%
    mutate(Comparison = case_when(
      cur_sample == 'ISG-6' & Visit %in% c('1', '2') ~ 'Group 1', # ISG-6
      cur_sample == 'ISG-6' & Visit == '5' ~ 'Group 2', # ISG-6
      cur_sample == 'ISG-34' & Visit == '1' ~ 'Group 1', # ISG-34
      cur_sample == 'ISG-34' & Visit %in% c('2', '3') ~ 'Group 2', # ISG-34
      TRUE ~ NA
    ))
  
  df_fc <- df %>% filter(!is.na(Comparison)) %>%
    group_by(Comparison) %>%
    summarise(across(all_of(columns), ~ mean(.x, na.rm = T))) %>%
    select(-Comparison) %>%
    t() %>%
    `colnames<-`(c('Group_1', 'Group_2')) %>%
    as.data.frame() %>%
    rownames_to_column(var = "name") %>%
    mutate(log2FC = log2(Group_2+1) - log2(Group_1+1)) %>%
    arrange(log2FC) %>%
    mutate(name = factor(name, levels=.$name))
  
  label_name <- rbind(df_fc %>% slice_max(log2FC, n=10),
                      df_fc %>% slice_min(log2FC, n=10))$name
  df_fc <- df_fc %>% mutate(label = if_else(name %in% label_name, name, NA))
  
  g <- ggplot(df_fc, aes(x=name, y = log2FC)) + geom_point() + 
    geom_label_repel(aes(label=label), box.padding = 0.5, max.overlaps = Inf) + 
    ylim(c(-10, 2)) +
    ggpubr::theme_pubclean() +
    theme(aspect.ratio=9/4, 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (cur_sample == 'ISG-34'){
    g <- g + ylab('log2FC(sample 2 and 3 / sample 1)') + ggtitle(cur_sample)
  } else if (cur_sample == 'ISG-6'){
    g <- g + ylab('log2FC(sample 5 / sample 1 and 2)') + ggtitle(cur_sample)
  }
}) %>%
  ggarrange(plotlist = ., ncol = 2) %>%
  ggexport(filename = 'figures/NFkB_foldchange.pdf')












