# trying out a new way to use tibbles
require(tidyverse)
require(broom)

# read in data nick sent me
nick <- read.csv("Downloads/Mean_SD_Final_2017.csv")
#make wavelength a factor
nick$wvl <- as.factor(nick$wvl)

# create a nest from grouping to support more general group-wise computation.
# Non-grouping variables are packaged into group-specific data frames held
# in a special variable called a list-column. We can then apply our computation
# to the components of the list (the data frames in the list) like an 'lapply' situation

nick_nested <- nick %>%
  group_by(lett, date, wvl) %>%
  nest()

# now write a function that we will subsequently apply to each data frame 
# we are going to first use the Kruskal-Wallis rank sum test 
# this is a non-parametric one-way anova on ranks. Sometimes called 
# an H-test. We can use it to compare  two or more independent samples of 
# equal or different sample sizes. The null hypothesis is that the mean ranks
# of the groups are the same.

kruskal_species_v_ref <- function(x) {
  kruskal.test(mean_ref ~ species, data = x)
}

#sanity check that our function is working
#kw <- kruskal.test(mean_ref~species, data = nick_nested[[1,"data"]])
#kwf <- kruskal_species_v_ref(nick_nested[[1,"data"]])
#kw$statistic == kwf$statistic

# Now, let's apply the function to all of the data frames 
# and store the results inside nick_nested in another column

# this is what it would look like for just the first 10 data frames
kws <- map(nick_nested$data[1:10], kruskal_species_v_ref)

# now lest do it to all and store it in the nested tibbles using mutate 
nick_nested <- nick_nested %>%
  mutate(kws = map(data, kruskal_species_v_ref))

# and let's do the same for our pairwise wilcox tests

wilcox_species_v_ref <- function(x){
  pairwise.wilcox.test(x$mean_ref, x$species, p.adj = "holm")
}

# and our same sanity check as above. Note that I'm querying the nests differently but it still works. 
# I'm pulling the 50000th dataframe out
# pw <- pairwise.wilcox.test(nick_nested$data[[50000]]$mean_ref, nick_nested$data[[50000]]$species)
# pwf <- wilcox_species_v_ref(nick_nested$data[[50000]])
# pw$p.value==pwf$p.value

nick_nested <- nick_nested %>%
  mutate(w.pairs = map(data, wilcox_species_v_ref))

# Now, let's also just check our asusmptions with a parametric approach 
# (some people argue that anovas are pretty robust anyhow). 
# this is always a nice sanity check that the results you are seeing are "real"
# For this, we are going to make a couple of steps. We are going to
# calculate the the anovas, and while we are at it the confidence intervals 
# and then add them all to our nested tibble


anova_lm_species_v_ref <- function(x) {
  anova(lm(mean_ref ~ species, data = x))
}

nick_nested <-nick_nested %>%
  mutate(avs = map(data,anova_lm_species_v_ref))

conf_avs <- function(x){
  confint(lm(mean_ref~species, data = x))
}

nick_nested <- nick_nested %>%
  mutate(cfs = map(data, conf_avs))

# So now we have a huge nested tibble with lots of statistics
# We can use the tidy() function to the model objects from lm(), etc.
# to create a data frame from our model results
# we will apply tidy() to the stats for each date-wvl-lett group using the 
# same map() inside mutate() strategy
# there is probably a way to do this with substaintially fewer lines of code, but I'm lazy
# like you should incoporate steps below into steps above

nick_nested <-nick_nested %>%
  mutate(tidy.kruskal = map(kws, tidy))

nick_nested <- nick_nested %>%
  mutate(tidy.wilcox = map(w.pairs, tidy))

nick_nested <- nick_nested %>%
  mutate(tidy.anovas = map(avs, tidy))

nick_nested <- nick_nested %>% 
  mutate(tidy.conf = map(cfs, tidy))

# now our last step is to get this all back into ordinary tibbles
# there is probably a better way to do this part, but I haven't figured it out

nick_kruskal <- nick_nested %>%
  select(lett, date, wvl, tidy.kruskal) %>%
  unnest(tidy.kruskal)

nick_wilcox <- nick_nested %>%
  select(lett, date, wvl, tidy.wilcox) %>%
  unnest(tidy.wilcox)

nick_anova <- nick_nested %>%
  select(lett, date, wvl, tidy.anovas) %>%
  unnest(tidy.anovas) %>%
  drop_na()

nick_confint <- nick_nested %>% 
  select(lett, date, wvl, tidy.conf) %>%
  unnest(tidy.conf) %>%
  filter(.rownames !="(Intercept)")

# let's combine the krskal wallis test results and the anova results 
# (since that makes sense for comparison)

colnames(nick_kruskal) <- c("lett", "date", "wvl", "h", "kw_p", "kw_df", "method")
colnames(nick_anova) <- c("lett", "date", "wvl", "term", "av_df", "sumsq", "meansq", "chi", "anova_p")
nick_k_a <- full_join(nick_kruskal,nick_anova)

# and let's see how well the anova and the kruskal perform

plot(nick_k_a$kw_p,nick_k_a$anova_p)
# phew, nice sanity check

# then for safety's sake, let's write these all as csv files 
# so we don't have to re-run everything

write.csv(nick_k_a, "/Users/ehestir/ncsu_laptop/workspace/kruskal_anova_stats.csv")
write.csv(nick_wilcox, "/Users/ehestir/ncsu_laptop/workspace/wilcox_pairs_stats.csv")
write.csv(nick_confint, "/Users/ehestir/ncsu_laptop/workspace/conf_int.csv")
