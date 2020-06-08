############################ NESP LCSP ANALYSIS ################################
# Rebekah Persad and Amie MacDonald
# June 2020
#
################################################################################

# load packages
library(tidyverse)
library(lubridate)
library(scales)

# load data
sparrow_data <- read_csv("nesp_lcsp_data.csv")

# select needed variables
sparrow_data <- sparrow_data %>% 
  select(Locality, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected,
         DayCollected, TimeObservationsStarted, DurationInHours, NumberOfObservers,
         ObservationCount, SpeciesCode)

# add zeros for days where surveys occurred, but no sparrows of one species were recorded
sparrow_data <- sparrow_data %>% 
  pivot_wider(names_from = SpeciesCode, values_from = ObservationCount,
              values_fill = list(ObservationCount = 0)) %>% 
  pivot_longer(cols = c("NESP", "LCSP"), names_to = "SpeciesCode", values_to = "ObservationCount") %>% 
  mutate(SurveyDate = paste0(MonthCollected, "-", DayCollected)) # create date column
  #mutate(SurveyDOY = yday(SurveyDate)) # create julian day column

# calculate relative abundance per year
sparrows_year <- sparrow_data %>% 
  group_by(SpeciesCode, YearCollected) %>% 
  summarize(MeanCount = mean(ObservationCount),
            SDCount = sd(ObservationCount)) %>% 
  ungroup()

# set custom theme for all plots
theme_ajm <- function() {
  theme_classic() %+replace%
    theme(axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=12, colour = "black"),
          axis.title.y = element_text(size=14, angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text = element_text(size=12),
          strip.text.x = element_text(size=12, face = "bold"),
          plot.title = element_text(size = 16, hjust = 0, vjust = 1),
          panel.border = element_rect(size =0.5, fill = "transparent"),
          plot.margin = margin(10, 10, 10, 15))
}

# plot relative abundance
png(filename = paste0("mean_count_year.png"),
    width=6, height=4, units="in", res=600)

ggplot(sparrows_year, aes(x=YearCollected, y=MeanCount, fill = SpeciesCode)) +
  geom_errorbar(aes(ymin=MeanCount-SDCount, ymax=MeanCount+SDCount), width=0, size=0.5, colour="black") +
  geom_line(linetype = "dashed") +
  geom_point(size=4, shape = 21) +
  scale_fill_manual(values = c("#440154FF", "#35B779FF"),
                    labels = c("LeConte's Sparrow", "Nelson's Sparrow")) +
  scale_x_continuous(breaks = seq(2009, 2016, 1)) +
  coord_cartesian(ylim = c(0, 33)) +
  ylab("Mean daily estimated total") +
  theme_ajm () +
  theme(legend.position = c(0.75, 0.85),
        legend.title = element_blank(),
        axis.title.x = element_blank())
  
dev.off()

# plot counts over season
sparrow_data$SurveyDate <- as.Date(sparrow_data$SurveyDate, format = "%m-%d")

png(filename = paste0("counts_season.png"),
    width=6, height=4, units="in", res=600)

ggplot(sparrow_data, aes(x=SurveyDate, y=ObservationCount, colour = SpeciesCode)) +
  geom_point(shape=16, alpha=0.5) +
  geom_smooth() +
  scale_colour_manual(values = c("#440154FF", "#35B779FF"),
                    labels = c("LeConte's Sparrow", "Nelson's Sparrow")) +
  scale_x_date(date_breaks="10 days", labels = date_format("%d %b")) +
  coord_cartesian(ylim = c(0, 25)) +
  xlab("Survey date") +
  ylab("Daily estimated total") +
  theme_ajm () +
  theme(legend.position = c(0.75, 0.85),
        legend.title = element_blank())

dev.off()

# plot by time of day
int1 <- interval(hms("00:00:00"), hms("04:29:00"))

sparrow_data <- sparrow_data %>% 
  mutate(TimeOfDay = case_when(TimeObservationsStarted >=  ~ "Evening",
                               TimeObservationsStarted >= "04:30:00" & "10:30:00" ~ "Morning",
                               TimeObservationsStarted >= "10:31:00" & "16:00:00" ~ "Afternoon",
                               TimeObservationsStarted >= "16:01:00" & "23:59:00" ~ "Evening"))
