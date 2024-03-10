#' """ Analysis for James et al. 2024 GCB
#'     @author: Ryan James
#'     date: 3/25/23"""

library(tidyverse)
library(sf)
library(stars)
library(gstat)
library(automap)
library(exactextractr)

# Generate habitat maps from monitoring data ----
# load shape files for each basin and make empty grid
# Rankin shape file
ran = st_read('gis/Rankin.shp') %>% st_transform(32617)

# make empty stars raster with 10 m cell size
ran_grid = st_bbox(ran) %>%
  st_as_stars(dx = 10) %>%
  st_set_crs(32617) %>%
  st_crop(ran)

# Johnson shape file
jon = st_read('gis/Johnson.shp') %>% st_transform(32617)

# make empty stars raster with 10 m cell size
jon_grid = st_bbox(jon) %>%
  st_as_stars(dx = 10) %>%
  st_set_crs(32617) %>%
  st_crop(jon)

# Whipray shape file
whp = st_read('gis/Whipray.shp') %>% st_transform(32617)

# make empty stars raster with 10 m cell size
whp_grid = st_bbox(whp) %>%
  st_as_stars(dx = 10) %>%
  st_set_crs(32617) %>%
  st_crop(whp)

# Create mangrove habitat map for each basin 
man = st_read('gis/FlBayMan500mBuf.shp')

# Rankin 
man_ran =  man |> st_intersection(ran) |> 
  st_rasterize(ran_grid) %>% 
  mutate(Acres = if_else(Acres > 0, 1, 0))
# Johnson
man_jon =  man |> st_intersection(jon) |>
  st_rasterize(jon_grid) %>% 
  mutate(Acres = if_else(Acres > 0, 1, 0))

# Whipray
man_whp =  man |> st_intersection(whp) |>
  st_rasterize(whp_grid) %>% 
  mutate(Acres = if_else(Acres > 0, 1, 0))

# combine into tibble
df_man = tibble(BASIN = c('RAN', 'JON', 'WHP'),
                man_r = c(list(man_ran), list(man_jon), list(man_whp)))

# Seagrass, algae, and total leaf area
df = st_read("FHAP_pc.csv",
             options=c("X_POSSIBLE_NAMES=LON","Y_POSSIBLE_NAMES=LAT")) |> 
  st_set_crs(4326) %>% st_transform(32617) |> 
  group_by(BASIN, YEAR) |> 
  nest() |> 
  mutate(sg = map(data, \(data) 
                  autofitVariogram(data$SG ~ 1, as(data, "Spatial"))),
         alg = map(data, \(data) 
                   autofitVariogram(data$ALG ~ 1, as(data, "Spatial"))),
         lab = map(data, \(data) 
                   autofitVariogram(data$LAB ~ 1, as(data, "Spatial"))),
         sg_r = case_when(
           BASIN == 'RAN' ~ map2(data,sg, \(data,sg)
                                 krige(SG~1, locations = data,
                                       newdata = ran_grid, model = sg$var_model)),
           BASIN == 'JON' ~ map2(data,sg, \(data,sg)
                                 krige(SG~1, locations = data,
                                       newdata = jon_grid, model = sg$var_model)),
           BASIN == 'WHP' ~ map2(data,sg, \(data,sg)
                                 krige(SG~1, locations = data,
                                       newdata = whp_grid, model = sg$var_model))),
         alg_r = case_when(
           BASIN == 'RAN' ~ map2(data,alg, \(data,alg)
                                 krige(ALG~1, locations = data,
                                       newdata = ran_grid, model = alg$var_model)),
           BASIN == 'JON' ~ map2(data,alg, \(data,alg)
                                 krige(ALG~1, locations = data,
                                       newdata = jon_grid, model = alg$var_model)),
           BASIN == 'WHP' ~ map2(data,alg, \(data,alg)
                                 krige(ALG~1, locations = data,
                                       newdata = whp_grid, model = alg$var_model))),
         lab_r = case_when(
           BASIN == 'RAN' ~ map2(data,lab, \(data,lab)
                                 krige(LAB~1, locations = data,
                                       newdata = ran_grid, model = lab$var_model)),
           BASIN == 'JON' ~ map2(data,lab, \(data,lab)
                                 krige(LAB~1, locations = data,
                                       newdata = jon_grid, model = lab$var_model)),
           BASIN == 'WHP' ~ map2(data,lab, \(data,lab)
                                 krige(LAB~1, locations = data,
                                       newdata = whp_grid, model = lab$var_model))))

# create dataset of rasters for each year
df_r = df |>
  select(YEAR, BASIN, sg_r, alg_r, lab_r) |> 
  left_join(df_man, by = 'BASIN')
           
# saveRDS(df_r, 'rasters.rds')
# E-scape ----
# Calculate IEI values ----
# load mixing model data
iso = read_csv('mm3spp.csv')

# random points within each basin
# set sample number 
n = 50

# set radius length
radius = 300

# filter raster to only 2019 to use for IEI values 
r19 = df_r |> 
  filter(YEAR == 2019)

# create buffer around random points for each basin
buf_ran = st_sample(ran, n) %>%
  st_as_sf %>%
  st_transform(32617)%>% # utm 17
  mutate(site = row_number(),
         BASIN = 'RAN')%>%
  st_buffer(dist = radius)

buf_jon = st_sample(jon, n) %>%
  st_as_sf %>%
  st_transform(32617)%>% # utm 17
  mutate(site = row_number(),
         BASIN = 'JON')%>%
  st_buffer( dist = radius)

buf_whp = st_sample(whp, n) %>%
  st_as_sf %>%
  st_transform(32617)%>% # utm 17
  mutate(site = row_number(),
         BASIN = 'WHP')%>%
  st_buffer( dist = radius)

# calculate f_habitat for each buffer 
ran_fhab = tibble(site = buf_ran$site,
                  BASIN = 'RAN') |> 
  mutate(f_sg = exact_extract(terra::rast(r19[r19$BASIN == 'RAN',]$sg_r[[1]]), buf_ran, 'mean', progress = F)[[1]],
         f_man = exact_extract(terra::rast(r19[r19$BASIN == 'RAN',]$man_r[[1]]), buf_ran, 'mean', progress = F)[[1]],
         f_epi = exact_extract(terra::rast(r19[r19$BASIN == 'RAN',]$lab_r[[1]]), buf_ran, 'mean', progress = F)[[1]],
         f_alg = exact_extract(terra::rast(r19[r19$BASIN == 'RAN',]$alg_r[[1]]), buf_ran, 'mean', progress = F)[[1]])


jon_fhab = tibble(site = buf_jon$site,
                  BASIN = 'JON') |> 
  mutate(f_sg = exact_extract(terra::rast(r19[r19$BASIN == 'JON',]$sg_r[[1]]), buf_jon, 'mean', progress = F)[[1]],
         f_man = exact_extract(terra::rast(r19[r19$BASIN == 'JON',]$man_r[[1]]), buf_jon, 'mean', progress = F)[[1]],
         f_epi = exact_extract(terra::rast(r19[r19$BASIN == 'JON',]$lab_r[[1]]), buf_jon, 'mean', progress = F)[[1]],
         f_alg = exact_extract(terra::rast(r19[r19$BASIN == 'JON',]$alg_r[[1]]), buf_jon, 'mean', progress = F)[[1]])


whp_fhab = tibble(site = buf_whp$site,
                  BASIN = 'WHP') |> 
  mutate(f_sg = exact_extract(terra::rast(r19[r19$BASIN == 'WHP',]$sg_r[[1]]), buf_whp, 'mean', progress = F)[[1]],
         f_man = exact_extract(terra::rast(r19[r19$BASIN == 'WHP',]$man_r[[1]]), buf_whp, 'mean', progress = F)[[1]],
         f_epi = exact_extract(terra::rast(r19[r19$BASIN == 'WHP',]$lab_r[[1]]), buf_whp, 'mean', progress = F)[[1]],
         f_alg = exact_extract(terra::rast(r19[r19$BASIN == 'WHP',]$alg_r[[1]]), buf_whp, 'mean', progress = F)[[1]])

# join together 
fhab = bind_rows(ran_fhab, jon_fhab, whp_fhab) |> 
  mutate(across(f_sg:f_alg, ~if_else(.x < 0, 0, .x)))

# calculate IEI for each species
IEI = fhab |> 
  slice(rep(1:n(), each = length(unique(iso$Species)))) |> 
  mutate(Species = rep(unique(iso$Species), times = nrow(fhab))) |> 
  left_join(iso, by = 'Species') |> 
  mutate(iei_sg = Seagrass/f_sg,
         iei_man = Mangrove/f_man,
         iei_epi = Epiphytes/f_epi,
         iei_alg = Algae/f_alg,
         iei_man = ifelse(iei_man == Inf, NA, iei_man)) |> 
  group_by(Species) |> 
  summarize(across(iei_sg:iei_alg, \(x) median(x, na.rm = T)))

# Make E-scape for each species, year, and basin----

# set radius length for size of landscape foraging units
radius = 300

# make a grid of landscape foraging units
ran_grd = st_make_grid(ran, cellsize = radius*2) %>% 
  st_sf(geometry = .) %>%
  filter(lengths(st_intersects(., ran)) > 0) |> 
  mutate(BASIN = 'RAN')

jon_grd = st_make_grid(jon, cellsize = radius*2) %>% 
  st_sf(geometry = .) %>%
  filter(lengths(st_intersects(., jon)) > 0) |> 
  mutate(BASIN = 'JON')

whp_grd = st_make_grid(whp, cellsize = radius*2) %>% 
  st_sf(geometry = .) %>%
  filter(lengths(st_intersects(., whp)) > 0) |> 
  mutate(BASIN = 'WHP')

grd = bind_rows(ran_grd, jon_grd, whp_grd) |> 
  group_by(BASIN) |> 
  nest(grid = geometry)

# calculate the f_hab for each basin and year
f_hab = df_r |> 
  left_join(grd, by = 'BASIN') |> 
  mutate(f_sg = map2(sg_r, grid, \(r, grid) 
                     exact_extract(terra::rast(r), grid, 'mean', progress = F)[[1]]),
         f_epi = map2(lab_r, grid, \(r, grid) 
                     exact_extract(terra::rast(r), grid, 'mean', progress = F)[[1]]),
         f_man = map2(man_r, grid, \(r, grid) 
                     exact_extract(terra::rast(r), grid, 'mean', progress = F)[[1]]),
         f_alg = map2(alg_r, grid, \(r, grid) 
                     exact_extract(terra::rast(r), grid, 'mean', progress = F)[[1]])) |> 
  dplyr::select(YEAR, BASIN, grid, f_sg:f_alg) |> 
  unnest(grid:f_alg) |> 
  ungroup()|> 
  mutate(across(f_sg:f_alg, ~if_else(.x < 0, 0, .x)))

# calculate HRI for all forgaing units for each species
HRI = f_hab |>  
  slice(rep(1:n(), each = length(unique(IEI$Species)))) |> 
  mutate(Species = rep(unique(IEI$Species), times = nrow(f_hab))) |> 
  left_join(IEI, by = 'Species') |> 
  mutate(HRI = f_sg*iei_sg + f_man*iei_man + 
           f_epi*iei_epi + f_alg*iei_alg) |> 
  st_as_sf()

# save output of all E-scapes 
# st_write(HRI, 'gis/FLBay_E-scapeAll.shp')

# plot the before and after figure here------

# Statistical Analyses ----
# ANOVA----
library(ggeffects)
library(emmeans)
library(multcomp)
library(multcompView)

df = st_read('gis/FLBay_E-scapeAll.shp') |> 
  drop_na(HRI) |> 
  mutate(do = case_when(
    YEAR <= 2015 ~ 'Pre Die-off',
    YEAR > 2015 ~ 'Post Die-off'),
    do = factor(do, levels = c('Pre Die-off','Post Die-off'))) |> 
 as_tibble()


m = aov(HRI ~ do*Basin*Species, data = df)

summary(m)

# check region comparisons
model_means = emmeans(object = m,
                      ~ do*Basin*Species)

# add letters to each mean
model_means_cld = cld(object = model_means,
                      adjust = "bonf",
                      Letters = letters,
                      sort = F,
                      alpha = 0.05) |> 
  mutate(code = str_replace_all(.group, ' ', ''))

# show output
model_means_cld

df = left_join(df, model_means_cld, by = c('Species', 'Basin', 'do'))

d = df |> group_by(Species, Basin, do) |> 
  mutate(v = max(HRI+ 0.1))

# plot with letters

ggplot(d, aes(Basin, HRI, fill = do, label = code))+
  geom_boxplot()+
  geom_text(aes(y = v), position = position_dodge(width = 0.75), size = 5)+
  scale_fill_manual(values = c('#4f7942', '#badaaf'))+
  labs(x = NULL)+
  facet_wrap(~Species)+
  theme_bw()+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 18, colour = "black"), 
        axis.text.x = element_text(size = 16, colour = "black"), 
        plot.title = element_text(size = 18, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20))

# break point ----
library(strucchange)

df = st_read('gis/FLBay_E-scapeAll.shp') |>
  as_tibble() |>
  group_by(Species, BASIN, YEAR) |>
  summarize(hri = mean(HRI, na.rm = T), 
            .groups = 'drop') |> 
  arrange(YEAR) |> 
  group_by(Species, BASIN) |> 
  nest() |> 
  mutate(bp = map(data, \(data) breakpoints(hri~1, data = data, h = 3)),
         fstat = map(data, \(data) sctest(Fstats(hri~1, data = data))),
         f = map_dbl(fstat, 'statistic'),
         p = map_dbl(fstat, 'p.value'),
         ci = map2(data, bp, \(data, bp) {
           if(!is.na(bp$breakpoints)) {
             confint(bp, data = data)[1]$confint
           } else {
             NA
           }})) |> 
  unnest_wider(ci, simplify = T, names_sep = '_') |> 
  mutate(change = map_dbl(data, \(data) data$YEAR[ci_1[,2]]),
         lci = map_dbl(data, \(data) data$YEAR[ci_1[,1]]),
         uci = map_dbl(data, \(data) data$YEAR[ci_1[,3]])) |> 
  dplyr::select(Species, BASIN, f, p, change, lci, uci)
         

# correlation ----         
df = st_read('gis/FLBay_E-scapeAll.shp') |>
  as_tibble() 

# algae 
cor.test(x = df$f_alg, y = df$HRI, method = 'pearson')

# epipyhtes
cor.test(x = df$f_epi, y = df$HRI)

#mangrove
cor.test(x = df$f_man, y = df$HRI)

#seagrass
cor.test(x = df$f_sg, y = df$HRI,  method = 'pearson')

