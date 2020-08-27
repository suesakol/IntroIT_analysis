#############################################################
### Analysis -- Introductory it construction      ###########
### Sakol Suethanapornkul and Sarut Supasiraprapa ###########
#############################################################

library(tidyverse)
library(stringr)
library(rsample)
library(ggrepel)

library(magrittr)

library(lme4)
library(lmerTest)


library(effects)
library(ggeffects)





set.seed(3478)










##### RQ 1: Associations between adjectives and variants #####




##### Preprocessing 

# Read file as tibble

corpusIntro <- read_csv("Data/Corpus_introIT.csv")


# Select columns for a distinctive collexeme analysis

corpusDat <- corpusIntro %>% 
  select(Variant, TaggedQueryItem)


# Prepare data by extracting adjectives from POS tag column
# (1) str_extract: in each line extract first instance of adjective (POS = JJ) 
# since we only have one _JJ per line, str_extract is sufficient
# plus, it vectorizes the matched patterns
# regex in "()" captures compound adj eg. far-fetched, well-known, and even
# words like O.K. with period between them by a group here --> [\w|.]
# (2) str_replace: substitute _JJ, _JJR, _JJS [regex = _JJ.*] with nothing

# NOTE: %>% is embedded inside mutate to outline steps taken inside this column

corpusDat <- corpusDat %>% 
  mutate(TaggedQueryItem = 
           str_extract(TaggedQueryItem, "([\\w|.]*)-?([\\w|.]*)_JJ.") %>% 
           str_replace(., "_JJ.*", "") %>% 
           str_to_lower(.)
         ) %>% 
  filter(!is.na(TaggedQueryItem))


# Convert superlative adjectives to stem & make spelling consistent
# Use case_when() for multiple if_else patterns and end with TRUE ~ existing vectors

corpusDat <- corpusDat %>% 
  mutate(TaggedQueryItem = case_when(TaggedQueryItem %in% c("best", "better") ~ "good",
                                     TaggedQueryItem == "cheaper" ~ "cheap",
                                     TaggedQueryItem == "clearer" ~ "clear",
                                     TaggedQueryItem == "costlier" ~ "costly",
                                     TaggedQueryItem %in% c("easier", "easiest") ~ "easy",
                                     TaggedQueryItem == "fairer" ~ "fair",
                                     TaggedQueryItem == "farfetched" ~ "far-fetched",
                                     TaggedQueryItem == "faster" ~ "fast",
                                     TaggedQueryItem %in% c("harder", "hardest") ~ "hard",
                                     TaggedQueryItem == "likelier" ~ "likely",
                                     TaggedQueryItem %in% c("o.k.", "ok") ~ "okay",
                                     TaggedQueryItem == "odder" ~ "odd",
                                     TaggedQueryItem == "pleasanter" ~ "pleasant",
                                     TaggedQueryItem %in% c("safer", "safest") ~ "safe",
                                     TaggedQueryItem %in% c("simpler", "simplest") ~ "simple",
                                     TaggedQueryItem == "sweeter" ~ "sweet",
                                     TaggedQueryItem == "tougher" ~ "tough",
                                     TaggedQueryItem == "truer" ~ "true",
                                     TaggedQueryItem %in% c("wiser", "wisest") ~ "wise",
                                     TaggedQueryItem == "worse" ~ "bad",
                                     TRUE ~ TaggedQueryItem
                                     )
         )




##### Checking frequencies

# Calculate frequency table of occurrences in Adj_that and Adj_to

corpusDat %>% 
  group_by(TaggedQueryItem) %>% 
  summarize(Adj_that = sum(Variant == "Adj_that"), Adj_to = sum(Variant == "Adj_to"))


# Inspect total number of items per variant

corpusDat %>% 
  count(Variant)


# Get type-token ratio of the complete data set

corpusDat %>% 
  group_by(Variant) %>% 
  summarize(type = n_distinct(TaggedQueryItem), 
            n = n(), 
            TTR = (type/n)*100
            )


# If words were distributed uniformly, how many instances per word in Adj-that?
# ANS = 58.61

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  group_by(Variant) %>% 
  summarize(uniform = sum(n)/length(TaggedQueryItem))


# Plot the distribution 
TopLable <- corpusDat %>% 
  filter(Variant == "Adj_that") %>%
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)
  

corpusDat %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = TaggedQueryItem), 
                           show.legend = FALSE, hjust = -0.2) +
  labs(x = "Rank order", y = "Frequencies") +
  geom_hline(yintercept = 58.61, color = "red", linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) +
  ggsave("Adj-that-fulldata.png", dpi = 300)

  
# Repeat the same step with the Adj-to variant
corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  group_by(Variant) %>% 
  summarize(uniform = sum(n)/length(TaggedQueryItem))


# Plot the distribution 

TopLable <- corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)


corpusDat %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = TaggedQueryItem), 
                           show.legend = FALSE, hjust = -0.2) +
  labs(x = "Rank order", y = "Frequencies") +
  geom_hline(yintercept = 46.37, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) +
  ggsave("Adj-to-fulldata.png", dpi = 300)
  
 



##### Drawing random samples from the full corpus

# Step 1: create 2 rows for Variants (one for Adj-that and  Adj-to) with group_by()
# Step 2: nest all the adjectives of each variant into the row with nest()
# At this stage, we get a tibble within a tibble (a list column)
# Step 3: ungroup and mutate the number of random samples we want to draw for each

# Step 4: [drawing random samples] with map2() and use sample_n() function
# map2(x,y,function, ...) operates on two inputs so we tell it to go to "data"
# a column where adjectives of each variant are and draw "n" samples from it
# with a function sample_n; we provide replace = FALSE to sample_n without replacement

# Step 5: drop the columns we don't need that is data and n and unnest the dataframe

corpusSam <- corpusDat %>% 
  group_by(Variant) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(n = c(12700, 12700)) %>% 
  mutate(samp = map2(data, n, sample_n, replace = FALSE)) %>% 
  select(-data, -n) %>%
  unnest(samp)


# Get type-token ratio of the balanced data set

corpusSam %>% 
  group_by(Variant) %>% 
  summarize(type = n_distinct(TaggedQueryItem), 
            n = n(), 
            TTR = (type/n)*100
            )


# Create a sample table for DCA with 'necessary' as an example
corpusSam %>% 
  mutate(Necessary = if_else(TaggedQueryItem == "necessary", 1, 0)) %>% 
  group_by(Variant) %>% 
  count(Necessary)

  



##### Checking to see if two data sets are correlated with kendall's tau 

# We join two tibbles which we create inside full_join
# 1st tibble comes from corpusDat (full dataset)
# 2nd tibble comes from corpusSam (balanced data set)
# sort = T in count() allows us to sort by descending frequency 
# min_rank() creates a rank order (rank 1, rank 2, rank 3, etc)
# min_rank(var) ranks lowest value to highest, min_rank(-var) from highest to lowest


# Adj-that variant

Adj_that <- full_join(corpusDat %>%
                        filter(Variant == "Adj_that") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>% 
                        mutate(rank = min_rank(-n)),
                      corpusSam %>%
                        filter(Variant == "Adj_that") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>%
                        mutate(rank = min_rank(-n)),
                      by = "TaggedQueryItem", 
                      suffix = c("_ori", "_sam")
                      )

# Adj-to variant

Adj_to <- full_join(corpusDat %>%
                        filter(Variant == "Adj_to") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>% 
                        mutate(rank = min_rank(-n)),
                      corpusSam %>%
                        filter(Variant == "Adj_to") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>%
                        mutate(rank = min_rank(-n)),
                      by = "TaggedQueryItem", 
                      suffix = c("_ori", "_sam")
                    )


# %$% from magrittr package "explodes" the tibble to work with base function witout $ 

Adj_that %$% 
  cor.test(rank_ori, rank_sam, method = "kendall")                     

Adj_to %$% 
  cor.test(rank_ori, rank_sam, method = "kendall")  





##### Work with balanced sample

# Calculate basic descriptive information
corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>% 
  mutate(percent = (n/sum(n)*100)) %>% 
  head(10)


# How many hapax?
corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n == 1) %>% 
  summarize(total = n_distinct(TaggedQueryItem))


# Adj-to variant

corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>% 
  mutate(percent = (n/sum(n)*100)) %>% 
  head(10)


corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n == 1) %>% 
  summarize(total = n_distinct(TaggedQueryItem))





##### Checking distribution of corpus data with relative entropy

# Relative entropy
# dispersion measure for categorical data (e.g., adjectives in each variant)
# approximates 1 as distributions become more even and 0 for uneven distributions

# Create a function to calculate relative entropy
hrel <- function(x){
  percentage <- x/sum(x)
  find_hrel <- -sum(percentage * log(percentage))/log(length(percentage))
  find_hrel
  }


# Calculate entropy for the whole dataset of adj-that & adj-to

corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>% 
  summarize(hrel = hrel(n))

corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>% 
  summarize(hrel = hrel(n))




##### Plot the distribution in balanced sample
# Adj-that

corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  group_by(Variant) %>% 
  summarize(uniform = sum(n)/length(TaggedQueryItem))


TopLable <- corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)


corpusSam %>% 
  filter(Variant == "Adj_that") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = TaggedQueryItem), 
                           show.legend = FALSE, hjust = -0.2) +
  geom_hline(yintercept = 51.2, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
  ) +
  ggsave("Adj-that-balanceddata.png", dpi = 300)




# Adj-to

corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(Variant, TaggedQueryItem) %>% 
  count() %>% 
  group_by(Variant) %>% 
  summarize(uniform = sum(n)/length(TaggedQueryItem))

TopLable <- corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 5)


corpusSam %>% 
  filter(Variant == "Adj_to") %>% 
  group_by(TaggedQueryItem) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(ID = row_number()) %>% 
  head(n = 100) %>% 
  ggplot(aes(x = ID, y = n)) +
  geom_point(color = "black") +
  geom_line(stat = "identity") +
  ggrepel::geom_text_repel(data = TopLable, aes(label = TaggedQueryItem), 
                           show.legend = FALSE, hjust = -0.2) +
  geom_hline(yintercept = 30.8, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
  ) +
  ggsave("Adj-to-balanceddata.png", dpi = 300)




##### DCA analysis with script from Gries (2015)
# Export a balanced corpus to .txt
corpusSam_dat <- as.data.frame(corpusSam)
write_delim(corpusSam_dat, path = "BalancedCorpus.txt", delim = "\t")








##### RQ 2: Matching between production and corpus data #####


##### Corpus data: DCA

# Read in DCA table; convert strings to lower characters and remove trailing space, if any

dca_dat <- read_delim("./Data/DCAresults.txt", delim = "\t",
                    col_types = cols(Words = col_character())
                    )

dca_dat <- dca_dat %>% 
  mutate(Words = str_to_lower(Words) %>% 
           str_replace(., " ", "")
         )


# Mutate new variable: grand-mean center and standardize coll.strength

dca_dat <- dca_dat %>% 
  mutate(Coll_strength_center = round((Coll_strength - mean(Coll_strength))/sd(Coll_strength), 
                                      digits = 3)
         )





##### Production data: Elicitation task

# Read in production data
# Convert strings to lower characters and remove trailing space
# Fix coding for Variant column to make it consistent with dca_dat
# Mutate new column: adding order of each subject by extracting number from Frame column

produc_dat <- read_delim("./Data/Product_introIT.txt", delim = "\t")

produc_dat <- produc_dat %>% 
  mutate(Response = str_to_lower(Response) %>% 
           str_replace(., " ", ""),
         Variant = case_when(Variant == "Adj-that" ~ "Adj_that",
                             Variant == "Adj-to" ~ "Adj_to")
         ) 


produc_dat <- produc_dat %>% 
  mutate(Trial_order = str_extract(Frame, pattern = "\\d.?")) %>% 
  relocate(Trial_order, .after = Frame) %>% 
  mutate(Trial_order = as.numeric(Trial_order)
         )


# Recode gender, NA to Not specified

produc_dat <- produc_dat %>% 
  mutate(Sex = if_else(is.na(Sex) == TRUE, "NotSpec", Sex))





##### Gather basic information

produc_dat %>% 
  group_by(Participant) %>% 
  slice(n = 1) %>% 
  group_by(L1) %>% 
  summarize(m = mean(Age, na.rm = TRUE), sd = sd(Age, na.rm = TRUE))


produc_dat %>% 
  group_by(Participant) %>% 
  slice(n = 1) %>% 
  group_by(L1, Sex) %>% 
  count()


produc_dat %>% 
  group_by(L1, past_participle, Misspell_wrong_incomplete) %>% 
  count()


# Inspect frequency counts (irrespective of variants)

produc_dat %>% 
  filter(past_participle == 0, Misspell_wrong_incomplete == 0) %>% 
  group_by(Response) %>% 
  count()


# Inspect number of types 
produc_dat %>% 
  filter(past_participle == 0, Misspell_wrong_incomplete == 0) %>% 
  select(Variant, Response) %>% 
  group_by(Variant) %>% 
  summarize(type = n_distinct(Response)) 

produc_dat %>% 
  filter(past_participle == 0, Misspell_wrong_incomplete == 0) %>% 
  select(Variant, Response) %>% 
  summarize(type_total = n_distinct(Response))


# Calculate percentage of non-target responses

produc_dat %>% 
  group_by(L1, past_participle, Misspell_wrong_incomplete) %>% 
  count() %>% 
  ungroup() %>% 
  select(L1, n) %>% 
  group_by(L1) %>% 
  mutate(percent = n/sum(n))


# How many responses per variant per participant after removing non-target responses?

produc_trim <- produc_dat %>% 
  filter(past_participle == 0 & Misspell_wrong_incomplete == 0) %>% 
  group_by(L1, Participant, Variant) %>% 
  count()


# Basic descriptive summary

produc_trim %>% 
  group_by(L1, Variant) %>% 
  summarize(m = mean(n),
            sd = sd(n)
            )


# Run a linear mixed effect to see any difference in terms of L1 and variant?

mod_num <- lmer(n ~ 1 + L1 * Variant + (1 | Participant), data = produc_trim,
                contrasts = list(Variant = contr.sum, L1 = contr.sum)
                )





##### Main statistical analysis


# Join tibbles (production data and dca data)

dat <- left_join(produc_dat %>% 
                   filter(past_participle == 0, Misspell_wrong_incomplete == 0) %>% 
                   select(!c(past_participle, Misspell_wrong_incomplete)), 
                 dca_dat %>% 
                   select(Words, starts_with("Freq_raw"), Preference, starts_with("Coll")), 
                 by = c("Response" = "Words")
                 )


# Inspect responses that did not show any preference as identified from DCA

dat %>% 
  filter(Preference == "no_preference")


# Code whether each response was a distinctive collexeme of a given variant
# Remove items that were not in DCA from analysis; from 914 down to 783 rows

dat <- dat %>% 
  filter(across(Freq_raw_that:Coll_strength_center, ~ !is.na(.))) %>% 
  mutate(Match = case_when(Variant == Preference ~ 1, 
                           Variant %in% c("Adj_that", "Adj_to") & Preference == "no_preference" ~ 1,
                           TRUE ~ 0)
         )


dat <- dat %>% 
  mutate(Degree = as.factor(Degree),
         L1 = as.factor(L1),
         Variant = as.factor(Variant)
         )


# After cases are dropped, how many are left, by L1 and variant

dat %>% 
  group_by(L1, Variant) %>% 
  count() 


# Add order in which words appear, after responses not in DCA dropped

dat <- dat %>% 
  group_by(Participant, Variant) %>% 
  mutate(Final_order = row_number()) %>% 
  relocate(Final_order, .after = Trial_order)


# Check if number of items differed by L1 and variant
# Step 1: Create a dataframe

num_dat <- dat %>% 
  group_by(L1, Participant, Variant) %>% 
  count()


num_dat %>% 
  group_by(L1, Variant) %>% 
  summarize(m = mean(n),
            sd = sd(n)
            )


# Step 2: 

model1 <- lmer(n ~ L1 * Variant + (1 | Participant), 
               contrasts = list(Variant = contr.sum, L1 = contr.sum),
               data = num_dat)


summary(model1)



# Center trial order [Trial_order]
dat <- dat %>% 
  group_by(Participant, Variant) %>% 
  mutate(Trial_center = (Trial_order - mean(Trial_order))/sd(Trial_order))



# Check frequencies of adjectives in production data
# Combine those hapax legomena into one category "other"

freq_prod <- dat %>% 
  group_by(Response) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(Response_recode = if_else(n == 1, "other", Response)) %>% 
  select(-n)


# Join two tibbles: Response_recode used for random intercept for item

dat <- left_join(dat, freq_prod, by = "Response")


dat <- dat %>% 
  relocate(Response_recode, .after = Response)


# Create one score for AWEQ
# We use rowwise() as a grouping option to calculate summary statistics per row/across rows 
# We mutate a mean of AWEQ with c_across() [a variant of across()] that works with rowwise()
# Then ungroup() to remove a grouping and go back to the original dataset

dat <- dat %>% 
  rowwise() %>% 
  mutate(AWEQ = mean(c_across(starts_with("AWEQ")))
         ) %>% 
  ungroup() %>% 
  mutate(AWEQ_center = round((AWEQ - mean(AWEQ))/sd(AWEQ), digits = 3)
         ) %>% 
  select(-AWEQ) %>% 
  relocate(AWEQ_center, .after = AWEQ_PUB)
  




##### Intercept-only model

# Begin with a plot: correct proportion by participants, not exactly aligned with bernouli trials
# This however shows there is substantial variability across subjects

dat %>% 
  group_by(Participant) %>% 
  summarize(n = n(),
            total = mean(Match),
            se = sd(Match)/sqrt(n)
            ) %>% 
  mutate(se = if_else(is.na(se), 0, se),
         ID = str_c("S",Participant)
         ) %>%
  ungroup() %>% 
  ggplot(aes(x = factor(ID, levels = unique(ID)), y = total)) +
  geom_pointrange(aes(ymin = total - se, ymax = total + se)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
        ) +
  ggsave("Proportion.png", width = 11, height = 11, dpi = 300)


# Run a model

m0 <- glmer(Match ~ 1 + (1 | Participant) + (1 | Response_recode),
            data = dat,
            family = binomial(link = "logit")
            )


# Calculate Intraclass correlation (ICC)
# formula var(intercept_adjust_participant)/ [var(intercept_adjust_participant) + [(pi)^2 + 3] ]

# Intercept-only model: SD of varying intercept for items and subjects

0.4292^2/((0.4292^2) + ((pi^2)/3))

0.1726^2/((0.1726^2) + ((pi^2)/3))

exp(0.55304)                      # Odds-ratio = 1.74
exp(0.55304)/(1+ exp(0.55304))    # Probability of a match = 63% across participants

# Since there is no predictor, deviance can be taken to suggest model misfit (Hox et al 2017)
# Deviance 1029.4





##### Model with level-1 predictor

m1 <- glmer(Match ~ 1 + Final_order + Variant + Coll_strength_center + 
              (1 | Participant) + (1 | Response_recode),
              contrasts = list(Variant = contr.sum),
              data = dat,
              family = binomial(link = "logit")
            )

# This model will warn us about a singular fit, caused by almost zero variance in the estimate of
# variance of participant-specific random intercept
# Another model is run with this term dropped, but no difference was found x^2(1) = 0, p < 1
# Random intercept coefficient is retained


# Add intra-level interaction and random slopes

m1.1 <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center + 
                (1 | Participant) + (1 | Response_recode),
              contrasts = list(Variant = contr.sum),
              data = dat,
              family = binomial(link = "logit")
              )


m1.2 <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center + 
                (1 + Variant * Coll_strength_center | Participant) + 
                (1 | Response_recode),
              contrasts = list(Variant = contr.sum),
              data = dat,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))
              )


# The model with correlated random slopes for the interaction b'w variant & coll.strength 
# and random intercepts for participant had a better fit than the one with no interaction
# X^2(4) = 17.394, p < 0.01
# Other models run (with one random slope, etc) not shown here





##### Model with level-2 predictor

m2 <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center + L1 + AWEQ_center + 
              (1 + Variant * Coll_strength_center| Participant) + 
              (1 | Response_recode),
            contrasts = list(Variant = contr.sum, L1 = contr.sum),
            data = dat,
            family = binomial(link = "logit"),
            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))
            )

# We did include graduate degree but the model with this predictor wasn't signficantly different from m2
# X^2(1) = 0, p < 1


# Add cross-level interaction
# From m1.2 to m2, variance of item-specific random intercepts did not chance
# Interaction terms added to the model to see if this variance is conditioned on level-2 predictor

m2.1 <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center * L1 + AWEQ_center + 
                (1 + Variant * Coll_strength_center| Participant) + 
                (1 | Response_recode),
              contrasts = list(Variant = contr.sum, L1 = contr.sum),
              data = dat,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))
              )

# Model (m2.1) had a better fit than the one with two-way interaction (m2) X^2(3) = 7.93, p < 0.05, 
# Deviance = 930.72. And a model with four-way interaction did not improve the fit X^2(7) = 2.61, p = 0.82


# Adjust slopes by items 

m2.2 <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center * L1 + AWEQ_center + 
                (1 + Variant * Coll_strength_center | Participant) + 
                (1 + L1| Response_recode),
              contrasts = list(Variant = contr.sum, L1 = contr.sum),
              data = dat,
              family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))
              )


# With random slopes on L1 or AWEQ by items, we observed (1) correlation between slopes & intercepts of -/+1
# This is a sign that our model is over-parameterized; variance of AWEQ was small (0.003)
# A model with item-specific adjustment on random slopes did not yield a better fit either, X^2(2) = 0.43, p = 0.8





##### Final Model

m_final <- glmer(Match ~ 1 + Final_order + Variant * Coll_strength_center * L1 + AWEQ_center + 
                   (1 + Variant * Coll_strength_center| Participant) + 
                   (1 | Response_recode),
                 contrasts = list(Variant = contr.sum, L1 = contr.sum),
                 data = dat,
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))
                 )


# m_final: deviance 930.72
# While dropping participant-specific random intercepts and slopes eliminated singular fits,
# a model without these terms was worse than the m_final X^2(10) = 34.29, p < 0.001

# When Coll.strength was dropped, m_final is better, X^2(7) = 17.692, p < 0.05
# When Variant was dropped, m_final is better, X^2(7) = 34.178, p < 0.001
# When uncorrelated, the two models did not differ X^2(1) < 1, p = 0.98 but m_final had smaller AIC


# The problem of this seems to stem from sd of by-participant random intercepts near 0
# Correlation of random slopes and intercepts is extremely unstable
# This is clearly shown with Bayesian model


contrasts(dat$Variant) <- contr.sum
contrasts(dat$L1) <- contr.sum

m_final_brm <- brms::brm(Match ~ 1 + Final_order + Variant * Coll_strength_center * L1 + AWEQ_center + 
                   (1 + Variant * Coll_strength_center| Participant) + 
                   (1 | Response_recode),
                 data = dat,
                 prior = c(prior(lkj(2), class = cor)),
                 family = bernoulli(link = "logit"),
                 warmup = 1000,
                 iter   = 4000,
                 chain  = 4,
                 inits  = 0,
                 seed   = 5690,
                 cores = parallel::detectCores(),
                 control = list(adapt_delta = 0.99)
                 )




##### Report the final model

# Extract fixed-effects estimates --- with fixef(mod_XXX) --- and 95% CI ---with confint(mod_XXX)
# And finally calculate odd ratio with an exponent -- exp --

sum_table <- confint(m_final, parm = "beta_", method = "Wald") %>% 
  as_tibble(rownames = "Parameters") %>% 
  mutate(Coeff = fixef(m_final)) %>% 
  relocate(Coeff, .after = Parameters) %>% 
  mutate(across(where(is.numeric), exp))


# Visualize the marginal effects (i.e., predicted probabilities of match)

# We will use ggpredict from ggeffects package to obtain fitted model
# For three-way interaction, we specify the terms in terms = c() using variable names from the model
# For continuous variables e.g., Coll_strength, we can get marginal effects at specified values

predict_df <- ggeffects::ggpredict(m_final, 
                                   terms = c("Coll_strength_center [-0.2,4.0,7.0, 10.0]",
                                             "Variant", 
                                             "L1"),
                                   )


# ggpredict return a tidy data structure, with predicted value calculated for the first term provided
# Here, we obtained predicted probabilities of match at various values of collostructional strength
# The other two terms are grouping variables

predict_df <- predict_df %>% 
  rename(coll_strength = x,
         variant = group,
         group = facet)


ggplot(data = predict_df, aes(x = coll_strength, y = predicted, group = variant)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  theme_bw() +  
  theme(text = element_text(size = 15),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  facet_grid(group ~ variant) 


ggsave("Interaction.png", dpi = 300)





##### RQ 2: Collexeme analysis comparing expert and novice data #####

# Randomly sample 914 words from balanced data to ensure the same number of items to production
# In production data, there were 440 responses for Adj_that and 474 for Adj_to

corpusSmall <- corpusSam %>% 
  group_by(Variant) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(n = c(440, 474)) %>% 
  mutate(samp = map2(data, n, sample_n, replace = TRUE)) %>% 
  select(-data, -n) %>%
  unnest(samp)


# Calculate type-token ratio of the smaller corpus

corpusSmall %>% 
  group_by(Variant) %>% 
  summarize(n = n(),
            type = n_distinct(TaggedQueryItem),
            ttr = (type/n)*100
            )


# Count overall types across the two variants

corpusSmall %>% 
  summarize(type = n_distinct(TaggedQueryItem)
            )


# Check if frequency ranks of the size-matched corpus and the balanced corpus correlated

Adj_that <- full_join(corpusSam %>% 
                        filter(Variant == "Adj_that") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>% 
                        mutate(rank = min_rank(-n)),
                      corpusSmall %>% 
                        filter(Variant == "Adj_that") %>% 
                        count(TaggedQueryItem, sort = TRUE) %>% 
                        mutate(rank = min_rank(-n)),
                      by = "TaggedQueryItem",
                      suffix = c("_balanced", "_small")
                      )

Adj_to <- full_join(corpusSam %>% 
                      filter(Variant == "Adj_to") %>% 
                      count(TaggedQueryItem, sort = TRUE) %>% 
                      mutate(rank = min_rank(-n)),
                    corpusSmall %>% 
                      filter(Variant == "Adj_to") %>% 
                      count(TaggedQueryItem, sort = TRUE) %>% 
                      mutate(rank = min_rank(-n)),
                    by = "TaggedQueryItem",
                    suffix = c("_balanced", "_small")
                    )


Adj_that %$% 
  cor.test(rank_balanced, rank_small, method = "kendall")

Adj_to %$% 
  cor.test(rank_balanced, rank_small, method = "kendall")



# DCA analysis on size-matched corpus and production corpus with script from Gries (2015)

# Size-matched corpus

corpusSmall_dat <- as.data.frame(corpusSmall)
write_delim(corpusSmall_dat, path = "SmallCorpus.txt", delim = "\t")


# Production corpus

producFinal_dat <- as.data.frame(produc_dat) %>% 
  filter(past_participle == 0, Misspell_wrong_incomplete == 0) %>% 
  select(Variant, Response) %>% 
  arrange(Variant)

write_delim(producFinal_dat, path = "Finalproduc.txt", delim = "\t")


# Read in files after conducting DCA
# Convert characters and remove trailing space and remove unwanted columns

dca_sizematch <- read_delim("./Data/DCA_SizeMatch.txt", delim = "\t",
                      col_types = cols(Words = col_character())
                      )

dca_sizematch <- dca_sizematch %>% 
  mutate(Words = str_to_lower(Words) %>% 
           str_replace(., " ", "")
         ) %>% 
  select(!c(matches(".*_exp_.*"), starts_with("delta"))
         )


dca_product <- read_delim("./Data/DCA_Production.txt", delim = "\t",
                        col_types = cols(Words = col_character())
                        )


dca_product <- dca_product %>% 
  mutate(Words = str_to_lower(Words) %>% 
           str_replace(., " ", "")
  ) %>% 
  select(!c(matches(".*_exp_.*"), starts_with("delta"))
         )


# Combine the two files together
# Then add another column that checks whether preference in corpus and production matched

dca_final <- full_join(dca_sizematch, dca_product, by = "Words", 
                       suffix = c(".corpus", ".produc")
                       ) %>% 
  mutate(match = if_else(Preference.corpus != Preference.produc, 0, 1)
         )


# Check how many unique adjective shared across the two data sets

dca_final %>% 
  filter(!is.na(match)) %>% 
  summarize(n = n_distinct(Words))

dca_final %>% 
  filter(match == 0) %>% 
  summarize(n = n_distinct(Words))

dca_final %>% 
  filter(match == 1) %>% 
  group_by(Preference.corpus) %>% 
  summarize(n = n_distinct(Words))


# Check ranks of distinctive collexemes, e.g., 

dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_to") %>% 
  select(Words, starts_with("Coll_strength")) %>% 
  arrange(-Coll_strength.produc)


# QQ plot suggests data are not normally distributed

dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_that") %>% 
  ggplot(aes(sample = Coll_strength.corpus)) + 
  geom_qq() +
  geom_qq_line() +
  theme_bw()


dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_to") %>% 
  ggplot(aes(sample = Coll_strength.corpus)) + 
  geom_qq() +
  geom_qq_line() +
  theme_bw() 


# Conduct correlation test
# Note: correlation between variables "rank.corpus/produc" is the same as correlation
#       between coll.strength.corpus/produc 

dca_final %>% 
  filter(match == 1) %>% 
  mutate(rank.corpus = min_rank(-Coll_strength.corpus),
         rank.produc = min_rank(-Coll_strength.produc)
         ) %$%
  cor.test(rank.corpus, rank.produc, method = "spearman")


dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_that") %>% 
  mutate(rank.corpus = min_rank(-Coll_strength.corpus),
         rank.produc = min_rank(-Coll_strength.produc)
         ) %$%
  cor.test(rank.corpus, rank.produc, method = "spearman")


dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_to") %>% 
  mutate(rank.corpus = min_rank(-Coll_strength.corpus),
         rank.produc = min_rank(-Coll_strength.produc)
         ) %$%
  cor.test(rank.corpus, rank.produc, method = "spearman")


# Plot the correlation

# Adj-that 
topitem <- dca_final %>% 
  filter(Words %in% c("clear", "likely", "true", "unlikely", "evident", "probable"))


dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_that") %>% 
  ggplot(aes(x = Coll_strength.corpus, Coll_strength.produc)) + 
  geom_point(size = 2, shape = 1, position = "jitter") +
  ggrepel::geom_text_repel(data = topitem,
                           aes(label = Words),
                           color = "black",
                           show.legend = FALSE,
                           nudge_x = -0.2,
                           nudge_y = -0.05) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, color = "black") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_blank(),
        panel.grid = element_blank()
        ) +
  ggsave("Cor_AdjThat.png", dpi = 300)


# Adj-to

topitem <- dca_final %>% 
  filter(Words %in% c("difficult", "necessary", "important", "easy", "crucial", "hard"))


dca_final %>% 
  filter(match == 1 & Preference.corpus == "Adj_to") %>% 
  ggplot(aes(x = Coll_strength.corpus, Coll_strength.produc)) + 
  geom_point(size = 2, shape = 1, position = "jitter") +
  ggrepel::geom_text_repel(data = topitem,
                           aes(label = Words),
                           color = "black",
                           show.legend = FALSE,
                           nudge_x = -0.2,
                           nudge_y = -0.05) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.2, color = "black") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_blank(),
        panel.grid = element_blank()
        ) +
  ggsave("Cor_AdjTo.png", dpi = 300)
