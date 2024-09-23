# run after data processing in mainv2.R

# longitudinal trend plot, test
source("report_functions.R")
library(ggpubr)
clinical_info_raw <- read_csv("data/ISG parameters combined table.csv")
# for (cur_sample in list(c("ISG-1", "ISG-1M", "ISG-1F"), "ISG-6", "ISG-34", "ISG-34F", "ISG-6M")) {
# cur_sample <- c("ISG-1", "ISG-1M", "ISG-1F")
cur_sample <- "ISG-6M"
healthy_range <- (ISG %>% filter(Group == "Healthy control"))$geomean %>% quantile(c(0.01, 0.99))
df <- ISG %>%
    prepare_timeline_table(
        cur_sample_sets = cur_sample,
        sample_info = sample_info_all %>% filter(`Patient ID` %in% cur_sample)
    )
if (!grepl("M|F", cur_sample[1])) {
    df <- df %>% mutate(connection_group = if_else(Group == "Patients" & !grepl("M|F", `Subject ID`), TRUE, FALSE))
} else {
    df <- df %>% mutate(connection_group = TRUE)
}

start_day <- df$`Date of sampling`[1]
df$timepoints <- difftime(df$`Date of sampling`, start_day, units = "days") %>% as.double()
# ISG-1
# vertical_line_label <- c(
#     "2023-03-09", "2023-04-06", "2023-05-08", "2023-06-05", "2023-07-04", "2023-07-17", "2023-08-16", "2023-09-11", "2023-10-10", "2023-11-07", "2023-12-07",
#     "2024-01-08", "2024-02-12", "2024-03-21", "2024-04-08", "2024-05-06"
# )
# ISG-6
# vertical_line_label <- c(
#     "2023-03-02", "2023-03-27", "2023-04-25", "2023-05-23", "2023-06-29", "2023-07-25", "2023-08-23", "2023-09-18",
#     "2023-10-16", "2023-11-13", "2023-12-11", "2024-01-08", "2024-02-05", "2024-03-22", "2024-04-16",
#     "2024-05-17", "2024-06-17"
# )
# ISG-34
# vertical_line_label <- c(
#     "2023-05-10", "2023-06-08", "2023-07-06", "2023-08-08", "2023-09-05", "2023-10-05", "2023-11-01",
#     "2023-12-04", "2024-01-09", "2024-02-08", "2024-03-06", "2024-04-04", "2024-05-06"
# )
# ISG-34F
# vertical_line_label <- c("2023-09-11", "2023-10-12", "2023-11-14", "2023-12-11", "2024-01-10", "2024-02-06")
# ISG-6M
vertical_line_label <- c(
    "2023-05-21", "2023-08-21", "2023-09-18", "2023-10-16", "2023-11-13", "2023-12-11",
    "2024-01-08", "2024-02-09"
)
# vertical_line_label <- c()

vertical_line <- difftime(as.Date(vertical_line_label), start_day, units = "days")

set.seed(42)
x_min <- min(df$timepoints, na.rm = TRUE) - 100
x_max <- max(df$timepoints, na.rm = TRUE) + 100
ggplot(df, aes(x = timepoints, y = geomean)) +
    geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
        fill = "lightblue", alpha = 0.05
    ) +
    geom_point() +
    geom_path(data = . %>% filter(connection_group), aes(group = connection_group), linetype = 2) +
    geom_vline(xintercept = vertical_line) +
    annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
    scale_x_continuous(breaks = unique(df$timepoints), name = "Days") +
    geom_text_repel(aes(label = label), size = 2.5) +
    # scale_y_continuous(trans="log1p") +
    annotate(
        geom = "text", x = x_min, y = healthy_range[2], hjust = 0, vjust = 1, label = "Range of healthy controls",
        color = "blue"
    ) +
    theme_bw() +
    theme(
        aspect.ratio = 3 / 15,
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank()
    )
ggsave(paste0("figures/", cur_sample[1], " IFN trendline days.pdf"))

vertical_line <- as.Date(vertical_line_label)
df <- df %>% mutate(`Date of sampling` = as.Date(`Date of sampling`))
x_min <- min(df$`Date of sampling`, na.rm = TRUE) - 100
x_max <- max(df$`Date of sampling`, na.rm = TRUE) + 100
ggplot(df, aes(x = `Date of sampling`, y = geomean)) +
    geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
        fill = "lightblue", alpha = 0.05
    ) +
    geom_point() +
    geom_path(data = . %>% filter(connection_group), aes(group = connection_group), linetype = 2) +
    geom_vline(xintercept = vertical_line) +
    annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
    geom_text_repel(aes(label = label), size = 2.5) +
    # scale_y_continuous(trans="log1p") +
    scale_x_date(breaks = unique(df$`Date of sampling`), date_labels = "%Y %b %d") +
    annotate(
        geom = "text", x = x_min, y = healthy_range[2], hjust = 0, vjust = 1, label = "Range of healthy controls",
        color = "blue"
    ) +
    theme_bw() +
    theme(
        aspect.ratio = 3 / 15,
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank()
    )
ggsave(paste0("figures/", cur_sample[1], " IFN trendline dates.pdf"))


# clinical information - trendline
clinical_info <- clinical_info_raw %>%
    filter(`Subject ID` == cur_sample[1]) %>%
    mutate(Datum = as.Date(Datum))
lapply(c("LKM", "CRP", "Kalpro", "Leukos", "Lymfos", "Neutros", "SpO2"), function(name) {
    clinical_info %>%
        filter(!is.na(!!sym(name))) %>%
        ggplot(aes_string(x = "Datum", y = name)) +
        geom_point() +
        geom_path(aes(group = `Subject ID`), linetype = 2) +
        geom_vline(xintercept = vertical_line) +
        annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
        scale_x_date(breaks = unique(clinical_info$Datum), date_labels = "%Y %b %d", limits = c(x_min, x_max)) +
        theme_bw() +
        theme(
            aspect.ratio = 3 / 15,
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.grid.minor.x = element_blank(),
            panel.background = element_blank()
        )
}) %>%
    ggarrange(plotlist = ., ncol = 1, nrow = 7, align = "h") %>%
    ggexport(filename = paste0("figures/", cur_sample[1], " parameters.pdf"), width = 10, height = 30)







# visualize

# df <- left_join(ISG %>% select(Sample, geomean),
#                 NFkb %>% select(Sample, geomean), by = "Sample") %>%
#   mutate(geomean = geomean.y/geomean.x)

# NFkb_ratio <- NFkb %>% select(-geomean) %>% left_join(df, by="Sample") %>%
#   filter(!(Group %in% c("Patients", "Patient IFN stimulation")) |
#            `Subject ID` == cur_sample) %>%
#   mutate(label = if_else(`Subject ID` == cur_sample,
#                          paste(`Subject ID`, Visit, sep="_"),
#                          NA))

# score_scatter_plot(NFkb_ratio) +
#   scale_y_continuous() +
#   ylab("NFkB to IFN score ratio")
# ggsave("figures/ISG-6_NFkB_ratio scatter.pdf")

# healthy_range <- (NFkb_ratio %>% filter(Group=="Healthy control"))$geomean %>% quantile(c(0.01,0.99))

# df <- NFkb_ratio %>%
#   prepare_timeline_table(
#     current_sample = cur_sample,
#     sample_info = sample_info_all %>% filter(`Patient ID` == cur_sample)
#   ) %>%
#   mutate(connection_group = if_else(Group == "Patients" & !grepl("M|F", Visit), TRUE, FALSE))

# df$timepoints <- difftime(df$`Date of sampling`, start_day, units="days") %>% as.double()
# vertical_line <- difftime(as.Date(vertical_line_label), start_day, units = "days")

# set.seed(42)
# x_min = min(df$timepoints, na.rm = T) -100
# x_max = max(df$timepoints, na.rm = T) +100
# ggplot(df, aes(x=timepoints, y=geomean)) +
#   geom_rect(aes(xmin = x_min, xmax = x_max, ymin = healthy_range[1], ymax = healthy_range[2]),
#             fill = "lightblue", alpha=0.05) +
#   geom_point() +
#   geom_path(data = . %>% filter(connection_group), aes(group=connection_group), linetype = 2) +
#   geom_vline(xintercept = vertical_line) +
#   annotate(geom = "text", x = vertical_line, y = max(df$geomean), label = vertical_line_label, angle = 90, vjust = 1, hjust = 1, size = 2) +
#   scale_x_continuous(breaks = unique(df$timepoints), name = "Days") +
#   geom_text_repel(aes(label=label), size = 2.5) +
#   # scale_y_continuous(trans="log1p") +
#   annotate(geom="text", x=x_min, y=healthy_range[2], hjust = 0, vjust = 1, label="Range of healthy controls",
#            color="blue")+
#   theme_bw() +
#   theme(aspect.ratio=3/15,
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.grid.minor.x = element_blank(),
#         panel.background = element_blank()) +
#   ylab("NFkB to IFN score ratio")

# ggsave("figures/ISG-34 NFkB_ratio trendline.pdf")

# differentially abundant ISGs
# cur_sample <- "ISG-6"
lapply(c("ISG-6", "ISG-34"), function(cur_sample) {
    columns <- na.omit(unlist(panel_info[, "NF-kB score"], use.names = F))
    df <- dat %>%
        select(all_of(columns), `Subject ID`, Visit) %>%
        filter(`Subject ID` == cur_sample) %>%
        mutate(Comparison = case_when(
            cur_sample == "ISG-6" & Visit %in% c("1", "2") ~ "Group 1", # ISG-6
            cur_sample == "ISG-6" & Visit == "5" ~ "Group 2", # ISG-6
            cur_sample == "ISG-34" & Visit == "1" ~ "Group 1", # ISG-34
            cur_sample == "ISG-34" & Visit %in% c("2", "3") ~ "Group 2", # ISG-34
            TRUE ~ NA
        ))

    df_fc <- df %>%
        filter(!is.na(Comparison)) %>%
        group_by(Comparison) %>%
        summarise(across(all_of(columns), ~ mean(.x, na.rm = T))) %>%
        select(-Comparison) %>%
        t() %>%
        `colnames<-`(c("Group_1", "Group_2")) %>%
        as.data.frame() %>%
        rownames_to_column(var = "name") %>%
        mutate(log2FC = log2(Group_2 + 1) - log2(Group_1 + 1)) %>%
        arrange(log2FC) %>%
        mutate(name = factor(name, levels = .$name))

    label_name <- rbind(
        df_fc %>% slice_max(log2FC, n = 10),
        df_fc %>% slice_min(log2FC, n = 10)
    )$name
    df_fc <- df_fc %>% mutate(label = if_else(name %in% label_name, name, NA))

    g <- ggplot(df_fc, aes(x = name, y = log2FC)) +
        geom_point() +
        geom_label_repel(aes(label = label), box.padding = 0.5, max.overlaps = Inf) +
        ylim(c(-10, 2)) +
        ggpubr::theme_pubclean() +
        theme(
            aspect.ratio = 9 / 4,
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    if (cur_sample == "ISG-34") {
        g <- g + ylab("log2FC(sample 2 and 3 / sample 1)") + ggtitle(cur_sample)
    } else if (cur_sample == "ISG-6") {
        g <- g + ylab("log2FC(sample 5 / sample 1 and 2)") + ggtitle(cur_sample)
    }
}) %>%
    ggarrange(plotlist = ., ncol = 2) %>%
    ggexport(filename = "figures/NFkB_foldchange.pdf")