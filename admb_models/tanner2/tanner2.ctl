#### SECTION B ####
#Fleet Key:
#1 = Directed Fishery (Tanner Crab)
#2 = Groundfish
#3 = Red King Crab
#4 = Snow Crab
#5 = NMFS Male
#6 = NMFS Female

#### B.1 Blocks ####

#Blocks correspond to changing estimated parameter value over time.
##EG. if you have different mortality in a cold phase vs. a warm phase
# Number of blocks: 1-6 for each fleet; 7 for natural mortality
  7
# Block details
#1965 1966  1967  1968  1969  1970  1971  1972  1973  1974  1975  1976  1977  1978  1979  1980  1981  1982  1983  1984  1985  1986  1987  1988  1989  1990  1991  1992  1993  1994  1995  1996  1997  1998  1999  2000  2001  2002  2003  2004  2005  2006  2007  2008  2009  2010  2011  2012  2013  2014  2015  2016  2017
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 4 5 6 1 1 1 1 1 1 1 1 7 8 9 10  11  1 1 1 12  12  12  12  12
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

#### B.2 Growth and Natural Mortality ####
  1                    # Type of maurity vector (1 = logistic, 2= spline)
  2                    # Basis for maturity ogive (1=age;2=length) ## CIA: added
  7                    # Type of growth curve (7= log-linear & gamma)
  0                    # 1 for offsets and 0 for raw
  1                    # 1= CV; 2= SD
  20                   # Number of size classes and animal can grow in 1 year
  1.000  1.000         # proportion to each platoon(m & f)
  0.000                # Proportion of within to between variance (platoons)

#Low= min value
#Hi= max value
#Init= starting value for parameter estimation
#Phase= estimation phase (negative value means do not estimate)
#Ptype= prior type (currently 0=no prior or 1= is prior)
#Pmean=prior mean
#PSD= prior standard deviation
#block= if you want multiple values over time (eg. diff mortality for warm phase and cold phase)
#UseDev= unique value in each year
#Dev 1 and Dev2= start year and end year for Dev
#SigDev= constrains Dev

# Low      Hi    Init    Phase Ptype PMean PSD Block UseDec  Dev1  2-Dec SigDev
##  Male
   0.00   0.35   0.349     -5     0   0.29  0.01 0 0 1900  2100  0.5  #1 natural mortality (male)
  -1.00   2.00   0.00     5     0   0.00  0.00 0 0 1901  2100  0.5  #2 natural mortality (male immature offset)
   2.00   182    37       -3    0   0     0    0 0 1902  2100  0.5 #3 length  1 (male)   ##CIA: work on growth parameters here; diff female and male
   2.00   182    177      -3    0   0     0    0 0 1903  2100  0.5 #4 length  2 (male)
   5.00   35.00  10       -1     0   0     0    0 0 1904  2100  0.5 #5 increment length  1 (male)
   0.00   5.00   5        -1     0   0     0    0 0 1905  2100  0.5 #6 increment length  2 (male)
   0.00   1.00   0.99      -5     0   0     0    0 0 1906  2100  0.5 #7 molt  prob  length  1 (male)
   0.00   1.00   0.005     -5     0   0     0    0 0 1907  2100  0.5 #8 molt  prob  length2 (male)
   0.10   0.75   0.5      -5     0   0     0    0 0 1908  2100  0.5 #9 CV  (male)
   0      1      0.00027  -1    0   0     0    0 0 1909  2100  0.5  #10 male  length-weight interce
   0      10     3.022134 -1    0   0     0    0 0 1910  2100  0.5  #11 male  length-weight slope
   0      400    120      -1    0   0     0    0 0 1911  2100  0.5  #12 Length-at-50% maturity
   0      0.4    0.088113 -1    0   0     0    0 0 1912  2100  0.5  #13 Maturity  slope
   -1     2      1        -1    0   0     0    0 0 1913  2100  0.5  #14 Maturity  asymptote
##  Female
   -1.00   1.00   0         5   0  -0.30 0.100 0 0 1915  2100  0.05 #15 natural mortality (female offs)
   -2.00   2.00   0         5   0  0.00  0.001 0 0 1916  2100  0.05 #16 natural mortality (female imature)
    2.00   182    27       -3   0  0     0     0 0 1917  2100  0.5 #17 length  1 (female)
    2.00   182    92       -3   0  0     0     0 0 1918  2100  0.5 #18 length  2 (female)
    1.00   20.00  9         -1   10 1     0     0 0 1904  2100  0.5 #19 increment length  1 (female)
    0.00   5.00   2         -1   2  1     0     0 0 1905  2100  0.5 #20 increment length  2 (female)
    0.00   1.00   0.9999    -6   1  0.9   0.1   0 0 1906  2100  0.5 #21 molt  prob  length  1 (female)
    0.00   0.01   0.005     -6   1  0.1   0.1   0 0 1907  2100  0.5 #22 molt  prob  length2 (female)
    0.10   0.75   0.5       -5   0  0     0     0 0 1908  2100  0.5 #23 CV  (female)
    0      1      0.000441  -1  0  0     0     0 0 1924  2100  0.5  #24 female  length-weight intercept
    0      10     2.898686  -1  0  0     0     0 0 1925  2100  0.5  #25 female  length-weight slope
    0      200    75        -1  0  0     0     0 0 1926  2100  0.5  #26 Length-at-50% maturity
    0      0.4    0.211052  -1  0  0     0     0 0 1927  2100  0.5  #27 Maturity  slope
    -1     2      1         -1  0  0     0     0 0 1928  2100  0.5  #28 Maturity  asymptote
    0      1      0.625     -1  0  0     0     0 0 1929  2100  0.5  #29 Phi for mature  biomasss
    -2     2      0         -1  0  0     0     0 0 1930  2100  0.5  #30 Sex ratio

#### B.3 Stock and Recruitment ####
# SR Parameters
#  min      max        init     phase
   1.50000  30.00000  19.00000   1     # R0
   -5.0000  5.00000    0        -2     # R1 offset ##CIA: Added, Andre's Dec. 2016 Update *
    .20000   1.00000   1.00000  -1     # Steepness
    .00000  10.00000   1.00000  -3     # SigmaR
    # !!!! =========================
    # !!!!  SRPars[5]: Cole added a prior to this in the .tpl file  !!!
  0.009000  100.51000  11.50000  6     # Alpha of rec fun *can set to earlier phase if needed
   1.0000   10.01000    4.00000  6     # Beta of rec fn   *can set to earlier phase if needed

# Last recruitment class
10

# Specifications for recruitment estimation
 1953                            # First Year
 2016                            # Last year
    1                            # Phase for recruiment
 1960                            # First Phase in
 1980                            # Flat top (part 1)
 2000                            # Flat top (part 2)
 2025                            # last Phase in
 1966                            # Last year of early devs
 0.1                             # variation of early devs

#### B.4 Fishing Mortality and Selectivity ####
# Specifications related to fishing mortality
     5                        # Number of years for tuning Fs
   5.0                        # Maximum fishing mortality rate

# Initial Fs
#  min      max        init     phase
    .00000   3.00000    .00000 -1                      # Fleet 1
    .00000   3.00000    .00000 -1                      # Fleet 2
    .00000   3.00000    .00000 -1                      # Fleet 3
    .00000   3.00000    .00000 -1                      # Fleet 4

#### B.5 Slectivity and Retention ####
# Selectivity parameter style (1= log-offset from females; 2= estimate parameters for males)
 2

# Selectvity
## type 0 = constant; 1 = logistic; 2 = double logistic; 3 = spline
## spline: if type 3, number of knots in spline function
## mirror can set one fleet's selectivity equal to another
## block= number of times to change param
## special = make fleet male/female
## retain= is it a retained catch
## retain block = prob changing over time
## surv disc = survival rate for discards

# Age selectvity specifications
# Type Spline-option Mirror, Block Special Retain Retain_block  survival_disc
     0             0      0      0       0      0            0         .679       # Fleet 1
     0             0      0      0       0      0            0         .679       # Fleet 1
     0             0      0      0       0      0            0         .20       # Fleet 2
     0             0      0      0       0      0            0         .20       # Fleet 2
     0             0      0      0       0      0            0         .679       # Fleet 3
     0             0      0      0       0      0            0         .679       # Fleet 3
     0             0      0      0       0      0            0         .679       # Fleet 4
     0             0      0      0       0      0            0         .679       # Fleet 4
     0             0      0      0       0      0            0         .00       # Fleet 5
     0             0      0      0       0      0            0         .00       # Fleet 5
     0             0      0      0       0      0            0         .00       # Fleet 6
     0             0      0      0       0      0            0         .00       # Fleet 6

# Length selectvity specifications
# Type Spline-option Mirror, Block Special Retain Retain_block
     1             0      0      1       0      1            0           # Fleet 1
     1             0      0      1       0      1            0           # Fleet 1
     1             0      0      2       0      0            0           # Fleet 2
     1             0      0      2       0      0            0           # Fleet 2
     1             0      0      3       0      0            0           # Fleet 3
     1             0      0      3       0      0            0           # Fleet 3
     2             0      0      4       0      0            0           # Fleet 4 (males)
     1             0      0      4       0      0            0           # Fleet 4 (females)
     4             0      0      5       1      0            0           # Fleet 5
     1             0      0      5       1      0            0           # Fleet 5
     1             0      0      6       2      0            0           # Fleet 6
     4             0      0      6       2      0            0            # Fleet 6

# Selectivity parameters
# Low Hi  Init  Phase Ptype Pmean PSD Block useDev  Dev1  Dev2  SigDev
# Selectivity - fleet 1
# fleet 1,  sex 1,  parameter 1 of  2
  120  160 150 2 0 135 10 1 0 1900  2100  0.005  # Male  Inflection  point ##CIA: work on this parameter
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
#fleet  1,  sex 1,  parameter 2 of  2
  0.01  1.001 0.79  2 0 0.5  0.2 1 0 1900  2100  0.005  # Male  Slope
  -2  2 0 3 # Male  Slope offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  #fleet  1,  sex 2,  parameter 1 of  2
  80  150 117 2 0 117 10 0 0 1900  2100  0.05  # Female  Inflection  point
  #fleet  1,  sex 2,  parameter 2 of  2
  0.1 0.4 0.14  2 0 0.14  0.1 0 0 1900  2100  0.05  # Female  Slope

# Selectivity - fleet 2
  40  250 140 2 0 140 20  2 0 1900  2100  0.05  # Male  Inflection  point
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  0.01  0.5 0.08  2 0 0.08  0.04  2 0 1900  2100  0.05  # Male  Slope
  -2  2 0 3 # Male  Slope offset
  -2  2 0 3 # Male  Slope offset
  40  120 70  2 0 70  10  2 0 1900  2100  0.05  # Female  Inflection  point
  -2  2 0 3 # Female  Inflection  point offset
  -2  2 0 3 # Female  Inflection  point offset
  0.005 0.5 0.03  2 0 0.03  0.02  2 0 1900  2100  0.05  # Female  Slope
  -2  2 0 3 # Female  Slope offset
  -2  2 0 3 # Female  Slope offset

# Selectivity - fleet 3
  95  200 150 2 0 150 20  3 0 1900  2100  0.05  # Male  Inflection  point
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  0.01  0.5 0.09  2 0 0.09  0.01  3 0 1900  2100  0.05  # Male  Slope
  -2  2 0 3 # Male  Slope offset
  -2  2 0 3 # Male  Slope offset
  50  170 125 2 0 125 30  3 0 1900  2100  0.05  # Female  Inflection  point
  -2  2 0 3 # Female  Inflection  point offset
  -2  2 0 3 # Female  Inflection  point offset
  0.05  0.5 0.2 2 0 0.2 0.04  3 0 1900  2100  0.05  # Female  Slope
  -2  2 0 3 # Female  Slope offset
  -2  2 0 3 # Female  Slope offset

# Low Hi  Init  Phase Ptype Pmean PSD Block useDev  Dev1  Dev2  SigDev
# Selectivity - fleet 4
  40.00 140.00 95.00  2 0 95  10  4 0 1900  2100  0.05  # Male  Peak
  -2  2 1 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2.00   2.00 -0.35  2 0 0  0.5 4 0 1900  2100  0.05  # 2nd parameter of  double-L
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2.00  15.00  6.00  2 0 5  5  4 0 1900  2100  0.05  # 3nd parameter of  double-L
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -2.00   10.00  9.00  2 0 5  5  4 0 1900  2100  0.05  # 4th parameter of  double-L
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
 -10.00  100.00 -0.50  2 0 0  0.5  4 0 1900  2100  0.05  # 5th parameter of  double-L
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset
  -10.0  100.00  0.50  2 0 0  0.5  4 0 1900  2100  0.05  # 6th parameter of  double-L
  -2  2 0 3 # Male  Inflection  point offset
  -2  2 0 3 # Male  Inflection  point offset

  50  170 125 2 0 125 30  3 0 1900  2100  0.05  # Female  Inflection  point
  -2  2 0 3 # Female  Inflection  point offset
  -2  2 0 3 # Female  Inflection  point offset
  0.05  0.5 0.2 2 0 0.2 0.1  3 0 1900  2100  0.05  # Female  Slope
  -2  2 0 3 # Female  Slope offset
  -2  2 0 3 # Female  Slope offset

# Selectivity - fleet 5 (NMFS male)
# Low Hi  Init  Phase Ptype Pmean PSD Block useDev  Dev1  Dev2  SigDev
  20  200 120  2  0 120 10  5 0 1900  2100  0.05 # Male  Inflection  point
  -2  2 0 3  # Male  Inflection  point offset
  -2  2 0 3  # Male  Inflection  point offset
  0.1 1.001 0.5 -2  0 0.5 0.25  5 0 1900  2100  0.05 # Male  Slope
  -2  2 0 3  # Male  Slope offset
  -2  2 0 3  # Male  Slope offset
  0.00  1 1.0 -1  0 0 0 0 0 1900  2100  0.005 # Male  asymptote
  20  200 80  -2  0 0 0 0 0 1900  2100  0.5 # Female  Inflection  point
  0.2 1.001 0.5 -2  0 0 0 0 0 1900  2100  0.5 # Female  Slope

# Selectivity - fleet 6 (NMFS female)
  20.0 200.000 100 -2  0  0    0  0 0 1900  2100  0.5   # Male  Inflection  point
   0.2   1.001 0.5 -2  0  0    0  0 0 1900  2100  0.5   # Male  Slope
  20.0 200.00 70.0  2  0 70   10  6 0 1900  2100  0.05  # Female  Inflection  point
  -2.0   2.000 0.0  3                                   # Female  Inflection  point offset
  -2.0   2.000 0.0  3                                   # Female  Inflection  point offset
   0.1   1.001 0.5  2  0 0.5 0.25 6 0 1900  2100  0.05  # Female  Slope
  -2.0   2.000 0.0  3                                   # Female  Slope offset
  -2.0   2.000 0.0  3                                   # Female  Slope offset
  0.00   1.000 1.0 -1  0  0   0   0 0 1900  2100  0.005 # Female asymptote

# Retention parameters
# Min      Max        Init       Phase
# Fleet 1
 10.00000 200.00000 140.00000     1                  #Male inflection point
   .00000   5.00000    0.5000     1                             #Male slope
   .00000   5.00000   1.00000    -1                        #Male proportion
 10.00000 100.00000  10.00000    -1                #Female inflection point
   .00000   5.00000    .10000    -1                           #Female slope
   .00000   5.00000    .00000    -1                      #Female proportion


#### B.6 Catchability and Additional Variance ####
# Priors for survey q (catchability)
## type 0 = catchability=1
## type 1 = estimate catchability, accounting for the prior
## type 2 = estimate catchability, ignoring prior
## mean of the prior for the log of catchability
## SD of the prior for log of catchability

# Type mean (of log)  SD  block
     0       .000  .0          1          # Survey 1
     0       .000  .0          1          # Survey 2
     0       .000  .0          1          # Survey 3
     0       .000  .0          1          # Survey 4
     1    -0.0555  .01         5          # Survey 5
     1    -0.0555  .01         6          # Survey 6

# Additional variance
     0       -5    5         -0.01     -2
     0       -5    5         -0.01     -2
     0       -5    5         -0.01     -2
     0       -5    5         -0.01     -2
     0       -5    5         -0.01     -2
     0       -5    5         -0.01     -2


# Sigma for catches
 .050

#### B.7 Weights ####
# Lambdas
##these are weighting factors that we will use later
#Fleet: 1    2     3      4     5     6
  1.000  1.000  1.000  1.000                    #Catch Likelihood
   .100  1.000  1.000  1.000                    #Discard likelihood
  1.000  1.000  1.000  1.000  5.000  5.000      #Index Likelihood
   .00    .00    .00    .00                     #Effort Likelihood
  #1.00   1.00   1.00   1.00                     #Effort Likelihood
   .100   .100   .100   .100   .500   .500      #Size likelihood
  1.000  1.000  1.000  1.000  1.000  1.000      #Age-size likelihood
   .000   .000   .000   .000   .000   .000      #Mean-size likelihood
  1.000                                         #Tagging data

#### B.8 Final Inputs ####
# Diagnostic level
 -2
# Stop of xx fucntion calls
 -1

# Checksum
12345678
