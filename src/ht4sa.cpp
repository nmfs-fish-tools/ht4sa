
#include <vector>


#define RCPP_NO_SUGAR
#include <Rcpp.h>

class parameter {
public:
    double value;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();
    bool is_random_effect = false;
    bool estimated = false;

    parameter(double value, double min,
            double max, bool estimated) :
    value(value), min(min), max(max),
    estimated(estimated) {

    }

    parameter(double value) {
        this->value = value;
    }

    parameter() {
        this->value = 0;
    }
};

class ht4sa_ss_recr_dist_pattern {
public:
    double GPattern;
    double month;
    double area;
    double age;
};

class ht4sa_ss_blocks_per_pattern {
public:
    double GPattern;
    double month;
    double area;
    double age;
};

class ht4sa_ss_MG_parms {
public:
    Rcpp::NumericVector LO;
    Rcpp::NumericVector HI;
    Rcpp::NumericVector INIT;
    Rcpp::NumericVector PRIOR;
    Rcpp::NumericVector PR_SD;
    Rcpp::NumericVector PR_type;
    Rcpp::NumericVector PHASE;
    Rcpp::NumericVector env_var_and_link;
    Rcpp::NumericVector dev_link;
};

class ht4sa_ss_control {
public:
    //warnings;
    double nseas;
    double N_areas;
    double Nages;
    double Nsexes;
    double Npopbins;
    double Nfleets;
    double Do_AgeKey;
    //fleetnames;
    //sourcefile;
    //type;
    //ReadVersion;
    //eof;
    double EmpiricalWAA;
    double N_GP;
    double N_platoon;
    double recr_dist_method;
    double recr_global_area;
    double recr_dist_read;
    double recr_dist_inx;
    //recr_dist_pattern;
    double N_Block_Designs;
    //blocks_per_pattern;
    //Block_Design;
    double time_vary_adjust_method;
    //time_vary_auto_generation;
    double natM_type;
    //natM;
    double GrowthModel;
    double Growth_Age_for_L1;
    double Growth_Age_for_L2;
    double Exp_Decay;
    double Growth_Placeholder;
    double N_natMparms;
    double SD_add_to_LAA;
    double CV_Growth_Pattern;
    double maturity_option;
    double First_Mature_Age;
    double fecundity_option;
    double hermaphroditism_option;
    double parameter_offset_approach;
    //MG_parms;
    //MGparm_seas_effects;
    double SR_function;
    double Use_steep_init_equi;
    double Sigma_R_FofCurvature;
    //SR_parms;
    //SR_parms_tv;
    double do_recdev;
    double MainRdevYrFirst;
    double MainRdevYrLast;
    double recdev_phase;
    double recdev_adv;
    double recdev_early_start;
    double recdev_early_phase;
    double Fcast_recr_phase;
    double lambda4Fcast_recr_like;
    double last_early_yr_nobias_adj;
    double first_yr_fullbias_adj;
    double last_yr_fullbias_adj;
    double first_recent_yr_nobias_adj;
    double max_bias_adj;
    double period_of_cycles_in_recr;
    double min_rec_dev;
    double max_rec_dev;
    double N_Read_recdevs;
    double F_ballpark;
    double F_ballpark_year;
    double F_Method;
    double maxF;
    double F_iter;
    //init_F;
    //Q_options;
    //Q_parms;
    //size_selex_types;
    //age_selex_types;
    //size_selex_parms;
    //age_selex_parms;
    //size_selex_parms_tv;
    double Use_2D_AR1_selectivity;
    double TG_custom;
    double DoVar_adjust;
    //Variance_adjustment_list;
    double maxlambdaphase;
    double sd_offset;
    //lambdas;
    double N_lambdas;
    double more_stddev_reporting;
};


//
//class ht4sa_ss_control {
//public:
//    double nseas;
//    double n_areas;
//    double n_ages;
//    double n_pop_bins;
//    double n_fleets;
//    double do_age_key;
//    Rcpp::StringVector fleet_names;
//    std::string source_file;
//    std::string type = "Stock_Synthesis_control_file";
//    bool eof = true;
//    double empirical_waa = 0;
//    double n_gp = 1.0;
//    double n_platoon = 1.0;
//    double recr_dist_method = 2.0;
//    double recr_global_area = 1.0;
//    double recr_dist_inx = 0.0;
//    ht4sa_ss_recr_dist_pattern recr_dist_patterns;
//    Rcpp::List time_vary_auto_generation;
//    double natM_type;
//    std::vector<Rcpp::DataFrame> natM;
//    double GrowthModel;
//    double Growth_Age_for_L1;
//    double Growth_Age_for_L2;
//    double Exp_Decay = -999;
//    double Growth_Placeholder;
//    double N_natMparms = 0;
//    double SD_add_to_LAA;
//    double CV_Growth_Pattern;
//    double maturity_option;
//    double First_Mature_Age;
//    double fecundity_option;
//    double hermaphroditism_option;
//    double parameter_offset_approach;
//
//    void add_recr_dist_pattern(ht4sa_ss_recr_dist_pattern pattern) {
//        this->recr_dist_patterns.push_back(pattern);
//    }
//
//
//};

RCPP_EXPOSED_CLASS(parameter)
RCPP_EXPOSED_CLASS(ht4sa_ss_control)

RCPP_MODULE(ht4sa) {
    Rcpp::class_<parameter>("parameter")
            .constructor()
            .constructor<double>()
            .constructor<parameter>()
            .field("value", &parameter::value)
            .field("min", &parameter::min)
            .field("max", &parameter::max)
            .field("is_random_effect", &parameter::is_random_effect)
            .field("estimated", &parameter::estimated);
    
    Rcpp::class_<ht4sa_ss_control>("ht4sa_ss_control")
            .constructor()
            //warnings;
            .field("nseas", &ht4sa_ss_control::nseas)
            .field("N_areas", &ht4sa_ss_control::N_areas)
            .field("Nages", &ht4sa_ss_control::Nages)
            .field("Nsexes", &ht4sa_ss_control::Nsexes)
            .field("Npopbins", &ht4sa_ss_control::Npopbins)
            .field("Nfleets", &ht4sa_ss_control::Nfleets)
            .field("Do_AgeKey", &ht4sa_ss_control::Do_AgeKey)
            //fleetnames;
            //sourcefile;
            //type;
            //ReadVersion;
            //eof;
            .field("EmpiricalWAA", &ht4sa_ss_control::EmpiricalWAA)
            .field("N_GP", &ht4sa_ss_control::N_GP)
            .field("N_platoon", &ht4sa_ss_control::N_platoon)
            .field("recr_dist_method", &ht4sa_ss_control::recr_dist_method)
            .field("recr_global_area", &ht4sa_ss_control::recr_global_area)
            .field("recr_dist_read", &ht4sa_ss_control::recr_dist_read)
            .field("recr_dist_inx", &ht4sa_ss_control::recr_dist_inx)
            //recr_dist_pattern;
            .field("N_Block_Designs", &ht4sa_ss_control::N_Block_Designs)
            //blocks_per_pattern;
            //Block_Design;
            .field("time_vary_adjust_method", &ht4sa_ss_control::time_vary_adjust_method)
            //time_vary_auto_generation;
            .field("natM_type", &ht4sa_ss_control::natM_type)
            //natM;
            .field("GrowthModel", &ht4sa_ss_control::GrowthModel)
            .field("Growth_Age_for_L1", &ht4sa_ss_control::Growth_Age_for_L1)
            .field("Growth_Age_for_L2", &ht4sa_ss_control::Growth_Age_for_L2)
            .field("Exp_Decay", &ht4sa_ss_control::Exp_Decay)
            .field("Growth_Placeholder", &ht4sa_ss_control::Growth_Placeholder)
            .field("N_natMparms", &ht4sa_ss_control::N_natMparms)
            .field("SD_add_to_LAA", &ht4sa_ss_control::SD_add_to_LAA)
            .field("CV_Growth_Pattern", &ht4sa_ss_control::CV_Growth_Pattern)
            .field("maturity_option", &ht4sa_ss_control::maturity_option)
            .field("First_Mature_Age", &ht4sa_ss_control::First_Mature_Age)
            .field("fecundity_option", &ht4sa_ss_control::fecundity_option)
            .field("hermaphroditism_option", &ht4sa_ss_control::hermaphroditism_option)
            .field("parameter_offset_approach", &ht4sa_ss_control::parameter_offset_approach)
            //MG_parms;
            //MGparm_seas_effects;
            .field("SR_function", &ht4sa_ss_control::SR_function)
            .field("Use_steep_init_equi", &ht4sa_ss_control::Use_steep_init_equi)
            .field("Sigma_R_FofCurvature", &ht4sa_ss_control::Sigma_R_FofCurvature)
            //SR_parms;
            //SR_parms_tv;
            .field("do_recdev", &ht4sa_ss_control::do_recdev)
            .field("MainRdevYrFirst", &ht4sa_ss_control::MainRdevYrFirst)
            .field("MainRdevYrLast", &ht4sa_ss_control::MainRdevYrLast)
            .field("recdev_phase", &ht4sa_ss_control::recdev_phase)
            .field("recdev_adv", &ht4sa_ss_control::recdev_adv)
            .field("recdev_early_start", &ht4sa_ss_control::recdev_early_start)
            .field("recdev_early_phase", &ht4sa_ss_control::recdev_early_phase)
            .field("Fcast_recr_phase", &ht4sa_ss_control::Fcast_recr_phase)
            .field("lambda4Fcast_recr_like", &ht4sa_ss_control::lambda4Fcast_recr_like)
            .field("last_early_yr_nobias_adj", &ht4sa_ss_control::last_early_yr_nobias_adj)
            .field("first_yr_fullbias_adj", &ht4sa_ss_control::first_yr_fullbias_adj)
            .field("last_yr_fullbias_adj", &ht4sa_ss_control::last_yr_fullbias_adj)
            .field("first_recent_yr_nobias_adj", &ht4sa_ss_control::first_recent_yr_nobias_adj)
            .field("max_bias_adj", &ht4sa_ss_control::max_bias_adj)
            .field("period_of_cycles_in_recr", &ht4sa_ss_control::period_of_cycles_in_recr)
            .field("min_rec_dev", &ht4sa_ss_control::min_rec_dev)
            .field("max_rec_dev", &ht4sa_ss_control::max_rec_dev)
            .field("N_Read_recdevs", &ht4sa_ss_control::N_Read_recdevs)
            .field("F_ballpark", &ht4sa_ss_control::F_ballpark)
            .field("F_ballpark_year", &ht4sa_ss_control::F_ballpark_year)
            .field("F_Method", &ht4sa_ss_control::F_Method)
            .field("maxF", &ht4sa_ss_control::maxF)
            .field("F_iter", &ht4sa_ss_control::F_iter)
            //init_F;
            //Q_options;
            //Q_parms;
            //size_selex_types;
            //age_selex_types;
            //size_selex_parms;
            //age_selex_parms;
            //size_selex_parms_tv;
            .field("Use_2D_AR1_selectivity", &ht4sa_ss_control::Use_2D_AR1_selectivity)
            .field("TG_custom", &ht4sa_ss_control::TG_custom)
            .field("DoVar_adjust", &ht4sa_ss_control::DoVar_adjust)
            //Variance_adjustment_list;
            .field("maxlambdaphase", &ht4sa_ss_control::maxlambdaphase)
            .field("sd_offset", &ht4sa_ss_control::sd_offset)
            //lambdas;
            .field("N_lambdas", &ht4sa_ss_control::N_lambdas)
            .field("more_stddev_reporting", &ht4sa_ss_control::more_stddev_reporting);
}

    //    Rcpp::class_<ht4sa_ss_recr_dist_pattern>("ht4sa_ss_recr_dist_pattern")
    //            .constructor()
    //            .field("GPattern", &ht4sa_ss_recr_dist_pattern::GPattern)
    //            .field("month", &ht4sa_ss_recr_dist_pattern::month)
    //            .field("area", &ht4sa_ss_recr_dist_pattern::area)
    //            .field("age", &ht4sa_ss_recr_dist_pattern::age);
    //}
