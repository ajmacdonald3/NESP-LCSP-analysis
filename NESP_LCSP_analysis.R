############################ NESP LCSP ANALYSIS ################################
# Rebekah Persad and Amie MacDonald
# June 2020
#
################################################################################

# load packages
library(tidyverse)
library(lubridate)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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
  mutate(SurveyDate = paste0(MonthCollected, "-", DayCollected)) %>%  # create date column
  mutate(SurveyDateYear = paste0(YearCollected, "-", MonthCollected, "-", DayCollected)) # create julian day column

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

# plot map of mean count
# set up basic things for maps 
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
# lakes <- ne_download(scale = "medium", type = 'lakes', category = 'physical',
# returnclass = "sf", destdir = "./map-data/lakes") # only need this first time downloading
lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
                 returnclass = "sf",
                 destdir = paste0(getwd(), "./map-data/lakes")) # use this if already downloaded shapefiles

# get mean count for each location
# Northbluff had multiple coordinates... use just one set for the site
sparrow_data <- sparrow_data %>% 
  mutate(DecimalLatitude = str_replace(DecimalLatitude, "51.4858193", "51.4879571")) %>% 
  mutate(DecimalLatitude = str_replace(DecimalLatitude, "51.7027", "51.4879571")) %>% 
  mutate(DecimalLongitude = str_replace(DecimalLongitude, "-80.4398775", "-80.4528165")) %>% 
  mutate(DecimalLongitude = str_replace(DecimalLongitude, "-80.567", "-80.4528165"))

# get mean and sample size for each location  
sparrows_loc <- sparrow_data %>%
  group_by(Locality, DecimalLatitude, DecimalLongitude, SpeciesCode) %>%
  summarize(SampleSize = length(unique(SurveyDateYear)),
            MeanCount = mean(ObservationCount)) %>% 
  filter(SampleSize > 10)

# map with point size based on mean count
# convert recv.locs to an sf object
sparrows_loc_sf <- st_as_sf(sparrows_loc, coords = c("DecimalLongitude", "DecimalLatitude"), crs = 4236)

sparrows_loc$DecimalLatitude <- as.numeric(sparrows_loc$DecimalLatitude)
sparrows_loc$DecimalLongitude <- as.numeric(sparrows_loc$DecimalLongitude)

# set map limits
xmin <- min(sparrows_loc$DecimalLongitude) - 1.25
xmax <- max(sparrows_loc$DecimalLongitude) + 0.75
ymin <- min(sparrows_loc$DecimalLatitude) - 0.5
ymax <- max(sparrows_loc$DecimalLatitude) + 0.75

# facet labels
labels <- c(LCSP = "LeConte's Sparrow", NESP = "Nelson's Sparrow")

# plot map
png(filename = "DET_map.png",
    width=8, height=5, units="in", res=600)

ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = sparrows_loc,
             aes(x = DecimalLongitude, y = DecimalLatitude,
                 colour = Locality, size = MeanCount), alpha = 0.7) +
  geom_text(data = sparrows_loc, aes(x = DecimalLongitude, y = DecimalLatitude, label = Locality),
            size = 3, hjust = "right", nudge_x = -0.15) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  scale_color_viridis_d(guide = "none") +
  scale_size_continuous(name = "Mean Daily Estimated Total", range = c(1,8)) +
  facet_wrap(. ~ SpeciesCode, labeller = labeller(SpeciesCode = labels)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 6),
        legend.position = "bottom", strip.background = element_rect(fill = "white"))

dev.off()

