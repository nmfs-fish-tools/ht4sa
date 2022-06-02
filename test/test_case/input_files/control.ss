#V3.30.18.00;_safe;_compile_date:_Sep 30 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis

#_data_and_control_files: BSH_n.dat // BSH_n.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
4 #_Nblock_Patterns
 1 2 3 1 #_blocks_per_pattern 
# begin and end years of blocks
 2006 2015
 2001 2005 2006 2015
 2011 2011 2012 2014 2015 2015
 1970 1970
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
 #_Age_natmort_by sex x growthpattern (nest GP in sex)
 0.785 0.488 0.37 0.306 0.267 0.24 0.221 0.207 0.196 0.187 0.18 0.175 0.171 0.167 0.164 0.161 0.159 0.157 0.156 0.155 0.154 0.153 0.152 0.151 0.151
 0.728 0.492 0.383 0.32 0.279 0.251 0.23 0.214 0.202 0.192 0.184 0.177 0.172 0.167 0.163 0.16 0.157 0.155 0.153 0.151 0.149 0.148 0.147 0.146 0.145
#
2 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
20 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
4 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
3 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
# Sex: 1  BioPattern: 1  Growth
 10 120 64.4 65 10 6 -4 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 40 410 244.6 400 10 6 -2 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.1 0.25 0.147 0.15 0.8 6 -4 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 -10 10 1 1 0.8 6 -4 0 0 0 0 0.5 0 0 # Richards_Fem_GP_1
 0.01 1 0.25 0.0834877 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 -3 3 -1.06443 0 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 5.388e-06 5.388e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Fem_GP_1
 -3 3.5 3.102 3.102 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 -3 300 156.6 55 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat50%_Fem_GP_1
 -3 3 -0.16 -0.16 0.8 6 -3 0 0 0 0 0.5 0 0 # Mat_slope_Fem_GP_1
 -3 50 45 45 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_scalar_Fem_GP_1
 -3 3 0 0 0.8 6 -3 0 0 0 0 0.5 0 0 # Eggs_exp_len_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
# Sex: 2  BioPattern: 1  Growth
 -3 3 0.059011 0 0.8 6 -3 0 0 0 0 0.5 0 0 # L_at_Amin_Mal_GP_1
 -3 3 0.068275 0 0.8 6 -2 0 0 0 0 0.5 0 0 # L_at_Amax_Mal_GP_1
 -3 3 -0.2 0 0.8 6 -3 0 0 0 0 0.5 0 0 # VonBert_K_Mal_GP_1
 -3 3 0 0 0.8 6 -3 0 0 0 0 0.5 0 0 # Richards_Mal_GP_1
 -3 3 0 0 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_young_Mal_GP_1
 -3 3 -1.4381 0 0.8 6 -3 0 0 0 0 0.5 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -3 3 3.293e-06 3.293e-06 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Mal_GP_1
 -3 3.5 3.225 3.225 0.8 6 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution  
 -4 4 0 0 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_GP_1
 -4 4 0 0 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_Area_1
 -4 4 4 0 99 0 -3 0 0 0 0 0.5 0 0 # RecrDist_month_1
#  Cohort growth dev base
 0.1 10 1 1 1 6 -1 0 0 0 0 0.5 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
7 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             3            20       10.2582             9            10             6          1          0          0          0          0          0          0          0 # SR_LN(R0)
          0.01             1         0.391           0.5           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_zfrac
          0.01            10             2             1           0.2             6         -4          0          0          0          0          0          0          0 # SR_surv_Beta
             0             2           0.3           0.6           0.8             6         -3          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             6         -1          0          0          0          0          0          4          1 # SR_regime
             0             0             0             0            99             0         -1          0          0          0          0          0          0          0 # SR_autocorr
# timevary SR parameters
 -5 5 -0.000919868 0 2.5 6 4 # SR_regime_BLK4add_1970
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1990 # first year of main recr_devs; early devs can preceed this era
2013 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -5 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 1 #_recdev_early_phase
 -1 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1978.99 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1992.32 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2012.46 #_last_yr_fullbias_adj_in_MPD
 2019.53 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.6094 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -10 #min rec_dev
 10 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1985E 1986E 1987E 1988E 1989E 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014F 2015F 2016F
#  -0.0317289 -0.0103961 0.305029 0.517253 -0.282686 0.177474 0.19891 -0.246105 0.151316 -0.177136 -0.000625463 0.1324 -0.0932803 -0.0588708 0.381204 0.241648 -0.0370074 -0.339964 0.162072 -0.299308 0.0424646 -0.351696 -0.303352 -0.295618 -0.0340382 0.111672 0.0979819 -0.444737 -0.31125 0 0 0
#
#Fishing Mortality info 
0.2 # F ballpark value in units of annual_F
2013 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
5 # max F (methods 2-4) or harvest fraction (method 1)
4  # N iterations for tuning in hybrid mode; recommend 3 (faster) to 5 (more precise if many fleets)
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 1
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 0.001 5 0.197224 0.01 99 6 1 # InitF_seas_1_flt_4F4_JPN_KK_SH
#
# F rates by fleet x season
# Yr:  1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# F1_MEX 0.00143158 0.00138598 0.00134193 0.00131275 0.00129458 0.00110281 0.0011607 0.0017206 0.00105459 0.00198445 0.00526439 0.00426397 0.00617658 0.00172089 0.00253752 0.0110234 0.0120965 0.0179705 0.00805176 0.0141744 0.0155338 0.0143694 0.0159118 0.00740628 0.00844231 0.0119888 0.0107179 0.0109923 0.00775924 0.00872534 0.00765519 0.0071276 0.00661835 0.011021 0.00822752 0.00852932 0.0106535 0.0146846 0.0157136 0.0162018 0.0130016 0.0134044 0.0144008 0.0180587 0.0180668 0.0180668
# F2_CAN 0 0 0 0 0 0 0 0 2.77688e-06 3.48869e-05 0 0 9.80352e-05 0 0.000290421 0.000467823 0.000870837 2.4772e-06 9.31122e-07 2.18181e-05 1.99214e-06 0 0 3.79161e-07 1.20604e-07 2.50007e-06 2.36317e-06 6.94472e-06 2.53951e-06 4.58891e-06 1.3819e-05 1.54187e-05 5.01756e-05 1.1193e-05 9.37351e-07 6.14173e-05 2.91977e-05 1.89838e-05 2.92743e-05 2.6574e-05 4.36999e-05 2.90406e-05 8.20338e-05 3.01964e-05 7.51963e-05 7.51963e-05
# F3_CHINA 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0014551 0.00132718 0.00111722 0.000985437 0.00119494 0.000707698 0.000845669 0.000494847 0.00114075 0.00144132 0.00260479 0.00327964 0.00250448 0.000995212 0.00241526 0.00241526
# F4_JPN_KK_SH 0.105059 0.0723051 0.0762915 0.0647158 0.068779 0.097708 0.136308 0.113345 0.126902 0.125672 0.114448 0.0736867 0.064691 0.0559714 0.0559279 0.0733699 0.06248 0.0580814 0.0510461 0.043381 0.0492341 0.0509619 0.0601296 0.0580265 0.0548835 0.052537 0.0623716 0.0585004 0.0659245 0.0856857 0.0912958 0.0720832 0.0670857 0.058035 0.0737436 0.0621123 0.0523554 0.0480305 0.0593876 0.0532342 0.0234619 0.0298888 0.0229816 0.0276375 0.0248913 0.0248913
# F5_JPN_KK_DP 0 0 0 0 0.00676047 0.0153871 0.0246267 0.0196871 0.019486 0.0407806 0.0616816 0.0685135 0.0832647 0.10608 0.0887099 0.103919 0.0588465 0.040105 0.0466444 0.0377178 0.0544698 0.0463429 0.0499722 0.0470427 0.025191 0.0368263 0.0216998 0.0144025 0.00983867 0.0059338 0.0044058 0.00431609 0.0100697 0.0299943 0.0102482 0.0172225 0.0329121 0.0104101 0.00428326 0.00223741 0.00130961 0.002781 0.0375557 0.0252193 0.000561897 0.000561897
# F6_JPN_ENY_SHL 0.00316471 0.0066703 0.00500056 0.00577954 0.00291734 0.0021537 0.00191934 0.00210152 0.00267314 0.00131011 0.00189838 0.00384887 0.0024929 0.0022122 0.00148192 0.00112399 0.00188032 0.00236042 0.00310059 0.00389709 0.00866325 0.00707711 0.0094609 0.0101451 0.0108552 0.0116403 0.0159148 0.0105072 0.0108398 0.0181623 0.0137559 0.0128237 0.0119971 0.0111119 0.0122877 0.0119021 0.00798664 0.0121078 0.00980314 0.00912602 0.0127605 0.0131164 0.00934342 0.00794719 0.0110486 0.0110486
# F7_JPN_ENY_DP 0 0 0 0 0.0214355 0.0385726 0.0520002 0.0517299 0.0786265 0.0885474 0.101874 0.0675726 0.0780906 0.0800613 0.0895305 0.0593318 0.0705019 0.117645 0.139325 0.0874376 0.07123 0.0563963 0.0712801 0.109243 0.138748 0.0767809 0.0826587 0.0754774 0.0497672 0.029803 0.0298421 0.0193454 0.0182433 0.0160238 0.0210825 0.0159655 0.016362 0.0156466 0.0123273 0.0382631 0.0565522 0.0165709 0.019347 0.0188685 0.0151108 0.0151108
# F8_JPN_LG_MESH 0 0 0.0408821 0.0407031 0.0407025 0.041069 0.0418177 0.0423935 0.0428322 0.0432526 0.0448197 0.0493704 0.0508073 0.0433853 0.036095 0.0408 0.040973 0.0356536 0.0318431 0.0346942 0.0347454 0.0306769 0.0060312 0.00487846 0.00403602 0.0038697 0.00465126 0.00473322 0.00643295 0.00515221 0.00479034 0.00505714 0.00963812 0.00859136 0.00989524 0.00899883 0.0103432 0.00774375 0.0102358 0.00787964 0.00361672 0.00414233 0.00402687 0.00382302 0.00354081 0.00354081
# F9_JPN_CST_Oth 0.00680734 0.0079591 0.00760911 0.0086443 0.00590735 0.0100329 0.00840154 0.0104945 0.0102692 0.00888214 0.00918829 0.00811014 0.00550353 0.0122342 0.0133098 0.0123907 0.0131046 0.0110825 0.00904098 0.0100161 0.00987437 0.00995185 0.00956043 0.0160594 0.0120212 0.00864805 0.00491416 0.00942616 0.00600232 0.0139056 0.00710583 0.00997792 0.011591 0.00881858 0.0188797 0.0180761 0.0219055 0.0208975 0.0190553 0.0156381 0.00441671 0.00614094 0.00641986 0.00411527 0.00491346 0.00491346
# F10_JPN_SM_MESH 0 0 0 0 0 0 0 0 0 0 0.40428 0.471453 0.507239 0.522685 0.534564 0.55452 0.51269 0.400438 0.595369 0.354371 0.30659 0.151 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# F11_IATTC 2.27752e-05 1.57498e-05 1.52492e-05 1.49176e-05 2.05955e-05 2.06409e-05 1.8042e-05 2.45361e-05 3.12009e-05 3.18021e-05 2.97423e-05 2.16629e-05 2.39402e-05 2.64752e-05 1.44177e-05 1.03604e-05 1.09719e-05 3.23113e-05 2.45032e-05 1.48423e-05 9.71776e-06 1.39735e-05 1.35728e-05 8.4258e-06 4.02015e-05 7.69252e-06 1.45426e-05 7.01487e-06 3.43177e-06 6.41805e-06 0 8.47179e-06 2.86881e-06 2.91485e-06 0 9.25423e-06 6.41003e-06 1.01157e-05 7.10542e-06 3.62537e-06 3.49599e-06 6.19969e-06 5.12711e-06 1.31289e-06 9.85104e-07 9.85104e-07
# F12_KOREA 0 0 3.81745e-08 1.73991e-07 1.68518e-05 0.000113496 0.000201509 6.47598e-05 0 0.000455192 1.16055e-06 0.00104215 0.000122871 0.000424763 0.000770062 0.000556778 0.00101298 0.000966519 0.00036373 0.000428177 0.00045147 0.000326544 0.000178043 0.000202021 0.000609584 0.00128302 0.00226482 0.00309138 0.00220166 0.00193249 0.000695515 0.00116758 0.00146009 0.00017302 0.000153088 7.55211e-05 0.000734096 0.000276378 0.000559335 0.00189602 0.00404843 0.00238502 0.00205498 0.00130259 0.000466177 0.000466177
# F13_NON_ISC 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000947793 0.000916934 0.0013654 0.00314511 0.00365692 0.00602514 0.00403504 0.0084586 0.00625395 0.0204128 0.0107247 0.0109752 0.0113851 0.00765064 0.00681163 0.00729362 0.0111587 0.0120247 0.00891886 0.00817694 0.00793272 0.00793272
# F14_USA_GIILL 0 0 0 0 0 0 0 7.5773e-06 1.92177e-05 4.68225e-05 0.00022419 0.000382922 0.000636917 0.000757683 0.000716993 0.00241584 0.00101984 0.000781517 0.000695953 0.00175179 0.000540395 0.000735025 0.000566529 0.000181443 0.000760948 0.000388287 0.000274345 0.000441058 0.00022457 0.000101522 8.08028e-05 5.41477e-05 7.39023e-05 4.91854e-05 1.57014e-05 1.01541e-05 9.51552e-05 5.1929e-05 1.95106e-05 1.18705e-05 1.13783e-05 1.76968e-05 7.01391e-06 7.20175e-06 3.577e-06 3.577e-06
# F15_USA_SPORT 0.000120501 0.000115895 0.000112098 0.000110021 0.000119842 0.000113138 0.000107937 0.000125025 0.000126837 0.000113154 0.000110057 6.8379e-05 0.000234385 0.000538797 0.00116285 0.00027628 0.00121441 0.00216324 0.000538276 0.000374966 0.000557641 0.000255898 0.000253589 0.000210866 0.000237796 0.00011877 0.000261485 4.62061e-05 8.31742e-05 0.000135362 5.83575e-05 2.25615e-05 5.4195e-05 1.96742e-05 1.57014e-05 1.35388e-05 1.76213e-05 1.11276e-05 1.17063e-05 1.18705e-05 3.79276e-06 7.43266e-06 4.55904e-06 5.7614e-06 6.7963e-06 6.7963e-06
# F16_USA_Lonline 0 0 0 0 0 0.000297925 0.000606569 0.000934226 0.00139827 0.00199158 0.00206149 0.00251501 0.00304193 0.00371842 0.00455369 0.00551488 0.00654421 0.00784604 0.00860092 0.0090728 0.00901028 0.0118483 0.0190884 0.0134028 0.00948604 0.0113102 0.0125059 0.0122542 0.0111926 0.00489619 0.00124474 0.000841015 0.000775115 0.000552061 0.000417795 0.000294109 0.000332922 0.00027727 0.000272171 0.000360525 0.000355958 0.000351852 0.000650162 0.000934608 0.0011045 0.0011045
# F17_TAIW_LG 1.37658e-05 1.36549e-05 2.69011e-06 0.000353226 0.000515494 2.0334e-05 0.000122169 0.000154377 3.66492e-05 0.000121151 0.00011129 1.43524e-05 1.47568e-05 3.05812e-06 0.000368681 0.000468563 0.000202227 3.86549e-05 0.00022491 0.000957891 0.00105789 0.000349724 0.000279473 5.54686e-05 0.00292096 0.0012325 0.00138474 0.00140232 0.00246936 0.00258347 0.00353185 0.00457318 0.00243693 0.00338336 0.00318818 0.00236604 0.00220375 0.00178379 0.00131696 0.00100612 0.00167444 0.00200223 0.0018453 0.00221982 0.00328456 0.00328456
# F18_TAIW_SM 0.0332168 0.0411039 0.0323459 0.0279629 0.0236921 0.0261241 0.0254442 0.0271267 0.0323194 0.0344472 0.030385 0.03444 0.0312284 0.0290766 0.0339732 0.0302528 0.0240976 0.0268883 0.0326661 0.0384442 0.043185 0.035159 0.0319484 0.025229 0.04332 0.0432138 0.0582155 0.047034 0.0534896 0.0744334 0.0317543 0.0344589 0.0296916 0.0326641 0.0472582 0.0355736 0.0351038 0.0371029 0.0443056 0.0392057 0.0446996 0.0534501 0.0233861 0.0353776 0.0245262 0.0245262
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
        19         1         0         0         0         1  #  S1_HW_DP
        21         1         0         0         0         1  #  S3_TAIW_LG
        23         1         0         0         0         1  #  S5_JPN_EARLY
        24         1         0         0         0         1  #  S6_JPN_LATE
        27         1         0         0         0         1  #  S9_SPC_OBS_TROPIC
        28         1         0         0         0         1  #  S10_MEX
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -9.20067             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S1_HW_DP(19)
           -25            25      -8.61177             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S3_TAIW_LG(21)
           -25            25      -8.52026             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S5_JPN_EARLY(23)
           -25            25       -8.7335             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S6_JPN_LATE(24)
           -25            25      -8.44813             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S9_SPC_OBS_TROPIC(27)
           -25            25      -9.37774             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_S10_MEX(28)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 4 0 # 1 F1_MEX
 5 0 0 1 # 2 F2_CAN
 24 0 4 0 # 3 F3_CHINA
 24 0 4 0 # 4 F4_JPN_KK_SH
 24 0 3 0 # 5 F5_JPN_KK_DP
 5 0 0 4 # 6 F6_JPN_ENY_SHL
 24 0 4 0 # 7 F7_JPN_ENY_DP
 24 0 4 0 # 8 F8_JPN_LG_MESH
 5 0 0 8 # 9 F9_JPN_CST_Oth
 24 0 0 0 # 10 F10_JPN_SM_MESH
 5 0 0 1 # 11 F11_IATTC
 5 0 0 3 # 12 F12_KOREA
 5 0 0 3 # 13 F13_NON_ISC
 24 0 4 0 # 14 F14_USA_GIILL
 5 0 0 14 # 15 F15_USA_SPORT
 24 0 4 0 # 16 F16_USA_Lonline
 24 0 4 0 # 17 F17_TAIW_LG
 5 0 0 17 # 18 F18_TAIW_SM
 5 0 0 16 # 19 S1_HW_DP
 5 0 0 16 # 20 S2_HW_SH
 5 0 0 17 # 21 S3_TAIW_LG
 5 0 0 18 # 22 S4_TAIW_SM
 5 0 0 4 # 23 S5_JPN_EARLY
 5 0 0 4 # 24 S6_JPN_LATE
 5 0 0 16 # 25 S7_JPN_RTV
 5 0 0 13 # 26 S8_SPC_OBS
 5 0 0 13 # 27 S9_SPC_OBS_TROPIC
 5 0 0 16 # 28 S10_MEX
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 11 0 0 0 # 1 F1_MEX
 11 0 0 0 # 2 F2_CAN
 11 0 0 0 # 3 F3_CHINA
 11 0 0 0 # 4 F4_JPN_KK_SH
 11 0 0 0 # 5 F5_JPN_KK_DP
 11 0 0 0 # 6 F6_JPN_ENY_SHL
 11 0 0 0 # 7 F7_JPN_ENY_DP
 11 0 0 0 # 8 F8_JPN_LG_MESH
 11 0 0 0 # 9 F9_JPN_CST_Oth
 11 0 0 0 # 10 F10_JPN_SM_MESH
 11 0 0 0 # 11 F11_IATTC
 11 0 0 0 # 12 F12_KOREA
 11 0 0 0 # 13 F13_NON_ISC
 11 0 0 0 # 14 F14_USA_GIILL
 11 0 0 0 # 15 F15_USA_SPORT
 11 0 0 0 # 16 F16_USA_Lonline
 11 0 0 0 # 17 F17_TAIW_LG
 11 0 0 0 # 18 F18_TAIW_SM
 11 0 0 0 # 19 S1_HW_DP
 11 0 0 0 # 20 S2_HW_SH
 11 0 0 0 # 21 S3_TAIW_LG
 11 0 0 0 # 22 S4_TAIW_SM
 11 0 0 0 # 23 S5_JPN_EARLY
 11 0 0 0 # 24 S6_JPN_LATE
 11 0 0 0 # 25 S7_JPN_RTV
 11 0 0 0 # 26 S8_SPC_OBS
 11 0 0 0 # 27 S9_SPC_OBS_TROPIC
 11 0 0 0 # 28 S10_MEX
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   F1_MEX LenSelex
            35           250       107.923            50             0             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F1_MEX(1)
           -15            15      -2.13702             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F1_MEX(1)
           -15            15       6.67147             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F1_MEX(1)
           -15            15       7.83583             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F1_MEX(1)
          -999          -999          -999             0             0             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F1_MEX(1)
          -999          -999          -999             0             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F1_MEX(1)
           -20           200      -14.6419           125            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F1_MEX(1)
           -15            15     -0.444108             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F1_MEX(1)
           -15            15     0.0154899             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F1_MEX(1)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F1_MEX(1)
           -15            15      0.608941             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F1_MEX(1)
# 2   F2_CAN LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F2_CAN(2)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F2_CAN(2)
# 3   F3_CHINA LenSelex
            35           250       176.483            50             0             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F3_CHINA(3)
           -15            15      -11.2974             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F3_CHINA(3)
           -15            15        6.5752             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F3_CHINA(3)
           -15            15       6.74023             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F3_CHINA(3)
          -999          -999          -999             0             0             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F3_CHINA(3)
          -999          -999          -999             0             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F3_CHINA(3)
           -20           200      -16.0333           125            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F3_CHINA(3)
           -15            15     -0.754683             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F3_CHINA(3)
           -15            15     -0.121267             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F3_CHINA(3)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F3_CHINA(3)
           -15            15      0.549526             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F3_CHINA(3)
# 4   F4_JPN_KK_SH LenSelex
            35           250       149.934            50             0             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F4_JPN_KK_SH(4)
           -15            15      -11.2542             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F4_JPN_KK_SH(4)
           -15            15       6.89704             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F4_JPN_KK_SH(4)
           -15            15       6.93663             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F4_JPN_KK_SH(4)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F4_JPN_KK_SH(4)
          -999          -999          -999             0             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F4_JPN_KK_SH(4)
           -20           200      0.674792             0            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F4_JPN_KK_SH(4)
           -15            15     -0.430054             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F4_JPN_KK_SH(4)
           -15            15     -0.953611             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F4_JPN_KK_SH(4)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F4_JPN_KK_SH(4)
           -15            15      0.362495             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F4_JPN_KK_SH(4)
# 5   F5_JPN_KK_DP LenSelex
            35           250       141.063            50             0             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F5_JPN_KK_DP(5)
           -15            15      -10.7303             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F5_JPN_KK_DP(5)
           -15            15       6.03132             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F5_JPN_KK_DP(5)
           -15            15         3.987             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F5_JPN_KK_DP(5)
          -999          -999          -999             0             0             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F5_JPN_KK_DP(5)
          -999          -999          -999             0             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F5_JPN_KK_DP(5)
           -80           200      -19.4448           125            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Peak_F5_JPN_KK_DP(5)
           -15            15     -0.395689             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Ascend_F5_JPN_KK_DP(5)
           -15            15       3.48731             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Descend_F5_JPN_KK_DP(5)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Male_Final_F5_JPN_KK_DP(5)
           -15            15      0.285421             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Male_Scale_F5_JPN_KK_DP(5)
# 6   F6_JPN_ENY_SHL LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F6_JPN_ENY_SHL(6)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F6_JPN_ENY_SHL(6)
# 7   F7_JPN_ENY_DP LenSelex
            35           250       167.883            50             0             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F7_JPN_ENY_DP(7)
           -15            15       -13.453             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F7_JPN_ENY_DP(7)
           -15            15       6.68519             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F7_JPN_ENY_DP(7)
           -15            15       6.62823             0             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F7_JPN_ENY_DP(7)
          -999          -999          -999             0             0             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F7_JPN_ENY_DP(7)
          -999          -999          -999             0             5             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F7_JPN_ENY_DP(7)
           -20           200      -10.7636           125            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Peak_F7_JPN_ENY_DP(7)
           -15            15     -0.458708             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F7_JPN_ENY_DP(7)
           -15            15     -0.312064             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F7_JPN_ENY_DP(7)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F7_JPN_ENY_DP(7)
           -15            15      0.589341             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F7_JPN_ENY_DP(7)
# 8   F8_JPN_LG_MESH LenSelex
            35           250         92.58           120             0             0          2          0          0          0          0        0.5          3          2  #  Size_DblN_peak_F8_JPN_LG_MESH(8)
           -15            15      -11.4556             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F8_JPN_LG_MESH(8)
           -15            15       6.99802             5             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F8_JPN_LG_MESH(8)
           -15            15       7.64907             5             0             0          4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F8_JPN_LG_MESH(8)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F8_JPN_LG_MESH(8)
          -999          -999          -999             0             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F8_JPN_LG_MESH(8)
           -20           200       199.934           125            50             0          4          0          0          0          0          0          3          2  #  SzSel_Fem_Peak_F8_JPN_LG_MESH(8)
           -15            15      0.424762             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F8_JPN_LG_MESH(8)
           -15            15     -0.427833             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F8_JPN_LG_MESH(8)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F8_JPN_LG_MESH(8)
           -15            15      0.394166             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F8_JPN_LG_MESH(8)
# 9   F9_JPN_CST_Oth LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F9_JPN_CST_Oth(9)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F9_JPN_CST_Oth(9)
# 10   F10_JPN_SM_MESH LenSelex
            35           250            50           120             0             0         -2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_F10_JPN_SM_MESH(10)
           -15            15            -9             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F10_JPN_SM_MESH(10)
           -15            15           5.5             5             0             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F10_JPN_SM_MESH(10)
           -15            15          6.85             5             0             0         -4          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F10_JPN_SM_MESH(10)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F10_JPN_SM_MESH(10)
          -999          -999          -999             0             4             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F10_JPN_SM_MESH(10)
# 11   F11_IATTC LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F11_IATTC(11)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F11_IATTC(11)
# 12   F12_KOREA LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F12_KOREA(12)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F12_KOREA(12)
# 13   F13_NON_ISC LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F13_NON_ISC(13)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F13_NON_ISC(13)
# 14   F14_USA_GIILL LenSelex
            28           250       71.7011            50             0             0          3          0          0          0          0        0.5          2          2  #  Size_DblN_peak_F14_USA_GIILL(14)
           -15            15       -1.7875             0             0             0          3          0          0          0          0        0.5          2          2  #  Size_DblN_top_logit_F14_USA_GIILL(14)
           -15            15       5.10292             0             0             0          3          0          0          0          0        0.5          2          2  #  Size_DblN_ascend_se_F14_USA_GIILL(14)
           -15            15       8.22599             0             0             0          3          0          0          0          0        0.5          2          2  #  Size_DblN_descend_se_F14_USA_GIILL(14)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F14_USA_GIILL(14)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F14_USA_GIILL(14)
           -20           200      -1.70798             0            50             0          4          0          0          0          0          0          2          2  #  SzSel_Fem_Peak_F14_USA_GIILL(14)
           -15            15      0.120667             4            50             0          4          0          0          0          0          0          2          2  #  SzSel_Fem_Ascend_F14_USA_GIILL(14)
           -15            15       -1.1752             4            50             0          4          0          0          0          0          0          2          2  #  SzSel_Fem_Descend_F14_USA_GIILL(14)
           -15            15             0             4            50             0          4          0          0          0          0          0          2          2  #  SzSel_Fem_Final_F14_USA_GIILL(14)
           -15            15      0.757923             4            50             0          5          0          0          0          0          0          2          2  #  SzSel_Fem_Scale_F14_USA_GIILL(14)
# 15   F15_USA_SPORT LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F15_USA_SPORT(15)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F15_USA_SPORT(15)
# 16   F16_USA_Lonline LenSelex
            35           250       175.593            50             0             0          3          0          0          0          0        0.5          1          2  #  Size_DblN_peak_F16_USA_Lonline(16)
           -15            15      -12.4935             0             0             0          3          0          0          0          0        0.5          1          2  #  Size_DblN_top_logit_F16_USA_Lonline(16)
           -15            15       7.13364             0             0             0          3          0          0          0          0        0.5          1          2  #  Size_DblN_ascend_se_F16_USA_Lonline(16)
           -15            15       6.90575             0             0             0          3          0          0          0          0        0.5          1          2  #  Size_DblN_descend_se_F16_USA_Lonline(16)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F16_USA_Lonline(16)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F16_USA_Lonline(16)
           -80           200      -25.1277             0            50             0          4          0          0          0          0          0          1          2  #  SzSel_Fem_Peak_F16_USA_Lonline(16)
           -15            15      -1.75007             4            50             0          4          0          0          0          0          0          1          2  #  SzSel_Fem_Ascend_F16_USA_Lonline(16)
           -15            15     -0.608215             4            50             0          4          0          0          0          0          0          1          2  #  SzSel_Fem_Descend_F16_USA_Lonline(16)
           -15            15             0             4            50             0          4          0          0          0          0          0          1          2  #  SzSel_Fem_Final_F16_USA_Lonline(16)
           -15            15      0.877276             4            50             0          5          0          0          0          0          0          1          2  #  SzSel_Fem_Scale_F16_USA_Lonline(16)
# 17   F17_TAIW_LG LenSelex
            35           250       214.908            50             0             0          2          0          1       2004       2014          6          0          0  #  Size_DblN_peak_F17_TAIW_LG(17)
           -15            15      -12.2676             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_F17_TAIW_LG(17)
           -15            15       7.26837             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_F17_TAIW_LG(17)
           -15            15       6.52454             0             0             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_F17_TAIW_LG(17)
          -999          -999          -999             0             0             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_F17_TAIW_LG(17)
          -999          -999          -999             0             5             0         -3          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_F17_TAIW_LG(17)
           -80           200      -13.4445             9            50             0          4          0          1       2004       2014          6          0          0  #  SzSel_Fem_Peak_F17_TAIW_LG(17)
           -15            15     0.0236907             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Ascend_F17_TAIW_LG(17)
           -15            15      0.811342             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Descend_F17_TAIW_LG(17)
           -15            15             0             4            50             0          4          0          0          0          0          0          0          0  #  SzSel_Fem_Final_F17_TAIW_LG(17)
           -15            15      0.482564             4            50             0          5          0          0          0          0          0          0          0  #  SzSel_Fem_Scale_F17_TAIW_LG(17)
# 18   F18_TAIW_SM LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_F18_TAIW_SM(18)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_F18_TAIW_SM(18)
# 19   S1_HW_DP LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S1_HW_DP(19)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S1_HW_DP(19)
# 20   S2_HW_SH LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S2_HW_SH(20)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S2_HW_SH(20)
# 21   S3_TAIW_LG LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S3_TAIW_LG(21)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S3_TAIW_LG(21)
# 22   S4_TAIW_SM LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S4_TAIW_SM(22)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S4_TAIW_SM(22)
# 23   S5_JPN_EARLY LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S5_JPN_EARLY(23)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S5_JPN_EARLY(23)
# 24   S6_JPN_LATE LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S6_JPN_LATE(24)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S6_JPN_LATE(24)
# 25   S7_JPN_RTV LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S7_JPN_RTV(25)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S7_JPN_RTV(25)
# 26   S8_SPC_OBS LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S8_SPC_OBS(26)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S8_SPC_OBS(26)
# 27   S9_SPC_OBS_TROPIC LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S9_SPC_OBS_TROPIC(27)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S9_SPC_OBS_TROPIC(27)
# 28   S10_MEX LenSelex
            -1           200             1            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_S10_MEX(28)
            -1           239            60            50            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_S10_MEX(28)
# 1   F1_MEX AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F1_MEX(1)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F1_MEX(1)
# 2   F2_CAN AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F2_CAN(2)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F2_CAN(2)
# 3   F3_CHINA AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F3_CHINA(3)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F3_CHINA(3)
# 4   F4_JPN_KK_SH AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F4_JPN_KK_SH(4)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F4_JPN_KK_SH(4)
# 5   F5_JPN_KK_DP AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F5_JPN_KK_DP(5)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F5_JPN_KK_DP(5)
# 6   F6_JPN_ENY_SHL AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F6_JPN_ENY_SHL(6)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F6_JPN_ENY_SHL(6)
# 7   F7_JPN_ENY_DP AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F7_JPN_ENY_DP(7)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F7_JPN_ENY_DP(7)
# 8   F8_JPN_LG_MESH AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F8_JPN_LG_MESH(8)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F8_JPN_LG_MESH(8)
# 9   F9_JPN_CST_Oth AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F9_JPN_CST_Oth(9)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F9_JPN_CST_Oth(9)
# 10   F10_JPN_SM_MESH AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F10_JPN_SM_MESH(10)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F10_JPN_SM_MESH(10)
# 11   F11_IATTC AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F11_IATTC(11)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F11_IATTC(11)
# 12   F12_KOREA AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F12_KOREA(12)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F12_KOREA(12)
# 13   F13_NON_ISC AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F13_NON_ISC(13)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F13_NON_ISC(13)
# 14   F14_USA_GIILL AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F14_USA_GIILL(14)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F14_USA_GIILL(14)
# 15   F15_USA_SPORT AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F15_USA_SPORT(15)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F15_USA_SPORT(15)
# 16   F16_USA_Lonline AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F16_USA_Lonline(16)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F16_USA_Lonline(16)
# 17   F17_TAIW_LG AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F17_TAIW_LG(17)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F17_TAIW_LG(17)
# 18   F18_TAIW_SM AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_F18_TAIW_SM(18)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_F18_TAIW_SM(18)
# 19   S1_HW_DP AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S1_HW_DP(19)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S1_HW_DP(19)
# 20   S2_HW_SH AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S2_HW_SH(20)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S2_HW_SH(20)
# 21   S3_TAIW_LG AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S3_TAIW_LG(21)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S3_TAIW_LG(21)
# 22   S4_TAIW_SM AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S4_TAIW_SM(22)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S4_TAIW_SM(22)
# 23   S5_JPN_EARLY AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S5_JPN_EARLY(23)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S5_JPN_EARLY(23)
# 24   S6_JPN_LATE AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S6_JPN_LATE(24)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S6_JPN_LATE(24)
# 25   S7_JPN_RTV AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S7_JPN_RTV(25)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S7_JPN_RTV(25)
# 26   S8_SPC_OBS AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S8_SPC_OBS(26)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S8_SPC_OBS(26)
# 27   S9_SPC_OBS_TROPIC AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S9_SPC_OBS_TROPIC(27)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S9_SPC_OBS_TROPIC(27)
# 28   S10_MEX AgeSelex
             0            40             0             1            99             0        -99          0          0          0          0        0.5          0          0  #  minage@sel=1_S10_MEX(28)
             1            40            24             3            99             0        -99          0          0          0          0        0.5          0          0  #  maxage@sel=1_S10_MEX(28)
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            35           250       125.107           120             0             0      4  # Size_DblN_peak_F8_JPN_LG_MESH(8)_BLK3repl_2011
            35           250       137.995           120             0             0      4  # Size_DblN_peak_F8_JPN_LG_MESH(8)_BLK3repl_2012
            35           250       152.427           120             0             0      4  # Size_DblN_peak_F8_JPN_LG_MESH(8)_BLK3repl_2015
           -80           200       5.42399           125            50             0      4  # SzSel_Fem_Peak_F8_JPN_LG_MESH(8)_BLK3repl_2011
           -80           200       15.1991           125            50             0      4  # SzSel_Fem_Peak_F8_JPN_LG_MESH(8)_BLK3repl_2012
           -80           200       8.38156           125            50             0      4  # SzSel_Fem_Peak_F8_JPN_LG_MESH(8)_BLK3repl_2015
            28           250        67.033            50             0             0      4  # Size_DblN_peak_F14_USA_GIILL(14)_BLK2repl_2001
            28           250       106.263            50             0             0      4  # Size_DblN_peak_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15      -2.01052             0             0             0      4  # Size_DblN_top_logit_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15      -2.45663             0             0             0      4  # Size_DblN_top_logit_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15       5.16316             0             0             0      4  # Size_DblN_ascend_se_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15       7.04233             0             0             0      4  # Size_DblN_ascend_se_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15       8.00452             0             0             0      4  # Size_DblN_descend_se_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15       8.20473             0             0             0      4  # Size_DblN_descend_se_F14_USA_GIILL(14)_BLK2repl_2006
           -20           200      -9.11408             0            50             0      4  # SzSel_Fem_Peak_F14_USA_GIILL(14)_BLK2repl_2001
           -20           200      -14.7349             0            50             0      4  # SzSel_Fem_Peak_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15     -0.674085             4            50             0      4  # SzSel_Fem_Ascend_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15     -0.218602             4            50             0      4  # SzSel_Fem_Ascend_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15     -0.612449             4            50             0      4  # SzSel_Fem_Descend_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15      -0.78566             4            50             0      4  # SzSel_Fem_Descend_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15             0             4            50             0      4  # SzSel_Fem_Final_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15             0             4            50             0      4  # SzSel_Fem_Final_F14_USA_GIILL(14)_BLK2repl_2006
           -15            15      0.674323             4            50             0      5  # SzSel_Fem_Scale_F14_USA_GIILL(14)_BLK2repl_2001
           -15            15      0.511123             4            50             0      5  # SzSel_Fem_Scale_F14_USA_GIILL(14)_BLK2repl_2006
            35           250       178.145            50             0             0      4  # Size_DblN_peak_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15      -11.1738             0             0             0      4  # Size_DblN_top_logit_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15       8.29167             0             0             0      4  # Size_DblN_ascend_se_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15       6.63012             0             0             0      4  # Size_DblN_descend_se_F16_USA_Lonline(16)_BLK1repl_2006
           -20           200      -9.11858             0            50             0      4  # SzSel_Fem_Peak_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15     -0.675171             4            50             0      4  # SzSel_Fem_Ascend_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15     -0.664103             4            50             0      4  # SzSel_Fem_Descend_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15             0             4            50             0      4  # SzSel_Fem_Final_F16_USA_Lonline(16)_BLK1repl_2006
           -15            15       1.00634             4            50             0      5  # SzSel_Fem_Scale_F16_USA_Lonline(16)_BLK1repl_2006
        0.0001             2           0.6           0.6           0.5             6      -5  # Size_DblN_peak_F17_TAIW_LG(17)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Size_DblN_peak_F17_TAIW_LG(17)_dev_autocorr
        0.0001             2           0.6           0.6           0.5             6      -5  # SzSel_Fem_Peak_F17_TAIW_LG(17)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # SzSel_Fem_Peak_F17_TAIW_LG(17)_dev_autocorr
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      2     5     1     4     1     0     0     0     0     0     0     0
#      5    60     2     3     2     0     0     0     0     0     0     0
#      5    66     5     3     2     0     0     0     0     0     0     0
#      5    85     8     2     2     0     0     0     0     0     0     0
#      5    86    10     2     2     0     0     0     0     0     0     0
#      5    87    12     2     2     0     0     0     0     0     0     0
#      5    88    14     2     2     0     0     0     0     0     0     0
#      5    91    16     2     2     0     0     0     0     0     0     0
#      5    92    18     2     2     0     0     0     0     0     0     0
#      5    93    20     2     2     0     0     0     0     0     0     0
#      5    94    22     2     2     0     0     0     0     0     0     0
#      5    95    24     2     2     0     0     0     0     0     0     0
#      5    98    26     1     2     0     0     0     0     0     0     0
#      5    99    27     1     2     0     0     0     0     0     0     0
#      5   100    28     1     2     0     0     0     0     0     0     0
#      5   101    29     1     2     0     0     0     0     0     0     0
#      5   104    30     1     2     0     0     0     0     0     0     0
#      5   105    31     1     2     0     0     0     0     0     0     0
#      5   106    32     1     2     0     0     0     0     0     0     0
#      5   107    33     1     2     0     0     0     0     0     0     0
#      5   108    34     1     2     0     0     0     0     0     0     0
#      5   109    35     0     0     0     0     1     1  2004  2014     6 -0.209646 0.305406 -0.0300147 -0.119742 -0.0289375 -0.10648 -0.200092 -0.0727251 -0.214978 -0.275437 -0.229577
#      5   115    37     0     0     0     0     2     1  2004  2014     6 -0.426445 1.95514 -1.09023 -0.0209612 0.0192866 0.295362 -0.409428 -0.896947 0.628128 -0.203256 0.121361
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      4      1      0.49
      4      5     0.408
      4     16      0.85
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 58 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 1 1 0 1
 1 2 1 0 1
 1 3 1 0 1
 1 4 1 0 1
 1 5 1 0 1
 1 6 1 0 1
 1 7 1 0 1
 1 8 1 0 1
 1 9 1 0 1
 1 10 1 0 1
 1 11 1 0 1
 1 12 1 0 1
 1 13 1 0 1
 1 14 1 0 1
 1 15 1 0 1
 1 16 1 0 1
 1 17 1 0 1
 1 18 1 0 1
 1 19 1 0 1
 1 20 1 0 1
 1 21 1 0 1
 1 22 1 0 1
 1 23 1 1 1
 1 24 1 1 1
 1 25 1 0 1
 1 26 1 0 1
 1 27 1 0 1
 1 28 1 0 1
 4 1 1 1 0
 4 2 1 0 0
 4 3 1 1 0
 4 4 1 1 0
 4 5 1 1 0
 4 6 1 0 0
 4 7 1 1 0
 4 8 1 1 0
 4 9 1 0 0
 4 10 1 0 0
 4 11 1 0 0
 4 12 1 0 0
 4 13 1 0 0
 4 14 1 1 0
 4 15 1 0 0
 4 16 1 1 0
 4 17 1 1 0
 4 18 1 0 0
 4 19 1 0 0
 4 20 1 0 0
 4 21 1 0 0
 4 22 1 0 0
 4 23 1 0 0
 4 24 1 0 0
 4 25 1 0 0
 4 26 1 0 0
 4 27 1 0 0
 4 28 1 0 0
 9 1 1 1 0
 12 1 1 1 0
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  0 #_CPUE/survey:_3
#  0 #_CPUE/survey:_4
#  0 #_CPUE/survey:_5
#  0 #_CPUE/survey:_6
#  0 #_CPUE/survey:_7
#  0 #_CPUE/survey:_8
#  0 #_CPUE/survey:_9
#  0 #_CPUE/survey:_10
#  0 #_CPUE/survey:_11
#  0 #_CPUE/survey:_12
#  0 #_CPUE/survey:_13
#  0 #_CPUE/survey:_14
#  0 #_CPUE/survey:_15
#  0 #_CPUE/survey:_16
#  0 #_CPUE/survey:_17
#  0 #_CPUE/survey:_18
#  0 #_CPUE/survey:_19
#  0 #_CPUE/survey:_20
#  0 #_CPUE/survey:_21
#  0 #_CPUE/survey:_22
#  1 #_CPUE/survey:_23
#  1 #_CPUE/survey:_24
#  0 #_CPUE/survey:_25
#  0 #_CPUE/survey:_26
#  0 #_CPUE/survey:_27
#  0 #_CPUE/survey:_28
#  1 #_lencomp:_1
#  0 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  0 #_lencomp:_6
#  1 #_lencomp:_7
#  1 #_lencomp:_8
#  0 #_lencomp:_9
#  0 #_lencomp:_10
#  0 #_lencomp:_11
#  0 #_lencomp:_12
#  0 #_lencomp:_13
#  1 #_lencomp:_14
#  0 #_lencomp:_15
#  1 #_lencomp:_16
#  1 #_lencomp:_17
#  0 #_lencomp:_18
#  0 #_lencomp:_19
#  0 #_lencomp:_20
#  0 #_lencomp:_21
#  0 #_lencomp:_22
#  0 #_lencomp:_23
#  0 #_lencomp:_24
#  0 #_lencomp:_25
#  0 #_lencomp:_26
#  0 #_lencomp:_27
#  0 #_lencomp:_28
#  1 #_init_equ_catch1
#  1 #_init_equ_catch2
#  1 #_init_equ_catch3
#  1 #_init_equ_catch4
#  1 #_init_equ_catch5
#  1 #_init_equ_catch6
#  1 #_init_equ_catch7
#  1 #_init_equ_catch8
#  1 #_init_equ_catch9
#  1 #_init_equ_catch10
#  1 #_init_equ_catch11
#  1 #_init_equ_catch12
#  1 #_init_equ_catch13
#  1 #_init_equ_catch14
#  1 #_init_equ_catch15
#  1 #_init_equ_catch16
#  1 #_init_equ_catch17
#  1 #_init_equ_catch18
#  1 #_init_equ_catch19
#  1 #_init_equ_catch20
#  1 #_init_equ_catch21
#  1 #_init_equ_catch22
#  1 #_init_equ_catch23
#  1 #_init_equ_catch24
#  1 #_init_equ_catch25
#  1 #_init_equ_catch26
#  1 #_init_equ_catch27
#  1 #_init_equ_catch28
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  1 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

