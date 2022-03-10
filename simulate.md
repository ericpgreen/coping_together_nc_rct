CT Trial Power Simulation
================

# Introduction

In this trial, families will be assigned to a treatment or an inactive
control arm. Families assigned to treatment will be nested in groups
that are assigned to facilitators. Data will be collected from all
eligible family members before (?) and after the intervention is
delivered to families assigned to the treatment arm.

![](diagram.png)<!-- -->

# Approach

We will conduct a Bayesian power simulation similar to what Kurz does in
[this
post](https://solomonkurz.netlify.app/post/bayesian-power-analysis-part-i/).

## Function to simulate data

The first step is to create a function to simulate data for this
structure. We’re using the
{[faux](https://debruine.github.io/faux/articles/sim_mixed.html)}
package.

``` r
# load the necessary packages
  library(faux)
  library(tidyverse)
  library(lme4)
  library(broom.mixed)
  library(brms)
  library(cmdstanr)
  library(tidybayes)
  library(sjPlot)
```

I’m unsure about:

1.  Setting group and facilitator columns to 0 for control arms
2.  Setting group/facilitator random effects to 0 for control arms
3.  Omitting random slopes
4.  Calculating dv: see `(b1 * treatment * post)`

``` r
#' Simulate data
#' @param seed  simulation seed
#' @param n_facilitator number of group facilitators
#' @param grp_per_fac_lo number of groups per facilitator, low end
#' @param grp_per_fac_hi number of groups per facilitator, high end
#' @param fam_per_gro_lo number of families per group, low end
#' @param fam_per_gro_hi number of families per group, high end
#' @param mem_per_fam_lo number of members per family, low end
#' @param mem_per_fam_hi number of members per family, high end
#' @param b0 intercept
#' @param b1 fixed effect of arm
#' @param u0l_sd random intercept SD for facilitators
#' @param u0g_sd random intercept SD for groups
#' @param u0f_sd random intercept SD for families
#' @param u0m_sd random intercept SD for members
#' @param u0t_sd random intercept SD for time

  sim <- function(n_facilitator, 
                  grp_per_fac_lo = 1,
                  grp_per_fac_hi = 3,
                  fam_per_gro_lo = 3,
                  fam_per_gro_hi = 4,
                  mem_per_fam_lo = 2,
                  mem_per_fam_hi = 5,
                  b0 = 0,
                  b1 = 0,       
                  u0l_sd = 0,   
                  u0g_sd = 0,   
                  u0f_sd = 0,   
                  u0m_sd = 0,   
                  u0t_sd = 0,   
                  sigma_sd = 0,
                  ... # helps the function work with pmap() 
                  ) {
    
# calculate nesting parameters
  n_groups_per_facilitator = sample(grp_per_fac_lo:grp_per_fac_hi, 
                                    n_facilitator, replace = TRUE)
  n_families_per_group = sample(fam_per_gro_lo:fam_per_gro_hi,
                                sum(n_groups_per_facilitator), replace = T)
  n_members_per_family = sample(mem_per_fam_lo:mem_per_fam_hi,
                                sum(n_families_per_group), replace = T)
  
# set up data structure
  df <- add_random(facilitator = n_facilitator) %>%
    add_random(group = n_groups_per_facilitator, 
               .nested_in = "facilitator") %>%
    add_random(family = n_families_per_group, .nested_in = "group") %>%
    add_random(member = n_members_per_family, .nested_in = "family") %>%
    add_within("member", time = c("pre", "post")) %>%
    add_between(.by = "group",
                arm = c("treatment", "control")) %>%
    add_recode("arm", "treatment", control = 0, treatment = 1) %>%
    add_recode("time", "post", pre = 0, post = 1) %>%
  # remove group/facilitator nesting for control arm
    mutate(facilitator = as.numeric(gsub("f", "", facilitator)),
           facilitator = case_when(
             arm == "control" ~ 0,
             TRUE ~ facilitator
    )) %>%
    mutate(group = as.numeric(gsub("g", "", group)),
           group = case_when(
             arm == "control" ~ 0,
             TRUE ~ group
    )) %>%
    mutate(family = as.numeric(gsub("f", "", family)),
           member = as.numeric(gsub("m", "", member))
           ) %>%
  # add random intercepts
    add_ranef("facilitator", u0l = u0l_sd) %>%
    add_ranef("group", u0g = u0g_sd) %>%
    add_ranef("family", u0f = u0f_sd) %>%
    add_ranef("member", u0m = u0m_sd) %>%
    add_ranef("time", u0t = u0t_sd) %>%
    add_ranef(sigma = sigma_sd) %>%
  # set group/facilitator random effects to 0 for control
    mutate(across(c(u0l, u0g), ~ case_when(arm == "control" ~ 0, 
                                           TRUE ~ .x))) %>%
  # TEMPORARY: Limit to 1 member per family
    distinct(family, time, .keep_all = TRUE) %>%
  # calculate DV
    mutate(dv = b0 + u0l + u0g + u0f + u0m + u0t + 
             (b1 * treatment * post) + sigma) %>%
  # reshape
    select(member, time, family, treatment, group, facilitator, dv) %>%
    pivot_wider(id_cols = c(member, family, treatment, group, facilitator),
                names_from = time,
                values_from = dv) %>%
    rename(y_pre = pre,
           y_post = post)
  }
```

## Simulate 1 time

Let’s imagine that the dv is a composite scale, specifically the mean of
10 items with possible values of 0 to 3. Thus, the dv can range from 0
to 3. Higher scores represent worse family functioning. Let’s also
imagine that we recruit distressed families, so the baseline mean is 2

I’m unsure about:

1.  How to set the random effects (I think family and member should
    remain 0 because we’re only dealing with 1 member per household at
    the moment)
2.  Does this look right?

``` r
  set.seed(8675309)
  data <- sim(# number of facilitators
                n_facilitator = 10, 
              # assume facilitators have 4 groups
                grp_per_fac_lo = 4, grp_per_fac_hi = 4,
              # assume groups have 3-4 families
                fam_per_gro_lo = 3, fam_per_gro_hi = 4,
              # assume families have 2-5 members
                mem_per_fam_lo = 2,
                mem_per_fam_hi = 5,
              # model parameters
                b0 = 2,             # grand mean
                b1 = -.3,           # treatment effect
                u0l_sd = 0,   
                u0g_sd = 0,   
                u0f_sd = 0,         # TEMP: set to 0 to look at 1 member/fam
                u0m_sd = 0,         # TEMP: set to 0 to look at 1 member/fam
                u0t_sd = 0,   
                sigma_sd = .7)
```

``` r
# look at the data
  data
```

    ## # A tibble: 91 × 7
    ##    member family treatment group facilitator y_pre y_post
    ##     <dbl>  <dbl>     <dbl> <dbl>       <dbl> <dbl>  <dbl>
    ##  1      1      1         1     1           1 1.65    1.82
    ##  2      4      2         1     1           1 2.16    2.02
    ##  3      7      3         1     1           1 3.06    1.67
    ##  4      9      4         0     0           0 1.67    2.59
    ##  5     12      5         0     0           0 1.35    2.61
    ##  6     16      6         0     0           0 2.09    1.10
    ##  7     18      7         1     3           1 1.12    2.90
    ##  8     22      8         1     3           1 0.789   1.39
    ##  9     27      9         1     3           1 2.08    1.56
    ## 10     32     10         1     3           1 1.41    1.23
    ## # … with 81 more rows

``` r
# families per arm
  data %>% count(treatment)
```

    ## # A tibble: 2 × 2
    ##   treatment     n
    ##       <dbl> <int>
    ## 1         0    43
    ## 2         1    48

Note: Because I let groups have varying numbers of families, there is
not a 1:1 allocation of families to arm.

``` r
  fit <- brm(y_post ~ treatment + y_pre + 
               (1 | group) + (1 | facilitator),
             data = data, 
             cores = parallel::detectCores(),
             backend = "cmdstanr")
```

``` r
  tab_model(fit)
```

<table style="border-collapse:collapse; border:none;">
<tr>
<th style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm;  text-align:left; ">
 
</th>
<th colspan="2" style="border-top: double; text-align:center; font-style:normal; font-weight:bold; padding:0.2cm; ">
y post
</th>
</tr>
<tr>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  text-align:left; ">
Predictors
</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">
Estimates
</td>
<td style=" text-align:center; border-bottom:1px solid; font-style:italic; font-weight:normal;  ">
CI (95%)
</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">
Intercept
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
2.16
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
1.51 – 2.83
</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">
treatment
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
-0.13
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
-0.77 – 0.46
</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; ">
y pre
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
-0.13
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:center;  ">
-0.32 – 0.05
</td>
</tr>
<tr>
<td colspan="3" style="font-weight:bold; text-align:left; padding-top:.8em;">
Random Effects
</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
σ<sup>2</sup>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
0.33
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
τ<sub>00</sub> <sub>facilitator</sub>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
0.02
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
τ<sub>00</sub> <sub>group</sub>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
0.03
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
ICC
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
0.13
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
N <sub>group</sub>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
14
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
N <sub>facilitator</sub>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
11
</td>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm; border-top:1px solid;">
Observations
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left; border-top:1px solid;" colspan="2">
91
</td>
</tr>
<tr>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; text-align:left; padding-top:0.1cm; padding-bottom:0.1cm;">
Marginal R<sup>2</sup> / Conditional R<sup>2</sup>
</td>
<td style=" padding:0.2cm; text-align:left; vertical-align:top; padding-top:0.1cm; padding-bottom:0.1cm; text-align:left;" colspan="2">
0.068 / 0.104
</td>
</tr>
</table>

# Next steps

1.  Priors!
2.  Attrition!
3.  Set up simulation to estimate power (including variants for effect
    size, random effects, families, etc)
4.  Talk about how to analyze data from multiple informants
