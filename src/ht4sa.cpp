
#include <vector>
#include <map>

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
    
    parameter(const parameter& other) :
    value(other.value), min(other.min), max(other.max), 
    is_random_effect(other.is_random_effect), 
    estimated(other.estimated) {
    }


    parameter(double value) {
        this->value = value;
    }
    
    

    parameter() {
        this->value = 0;
    }
};

enum ModelCategory {
    RECRUITMENT = 0,
    SELECTIVITY,
    GROWTH
};

class model_base {
public:
    static std::vector<model_base* > model_objects;
};
std::vector<model_base* > model_base::model_objects;

class recruitment_base : public model_base {
public:
    static std::map<uint32_t, recruitment_base*> recruitment_objects;
    static uint32_t id_g;
    uint32_t id;

    recruitment_base() {
        this->id = recruitment_base::id_g++;
        model_base::model_objects.push_back(this);
    }

};
std::map<uint32_t, recruitment_base*> recruitment_base::recruitment_objects;
uint32_t recruitment_base::id_g = 1;

class ricker : public recruitment_base {
public:
    uint32_t category = RECRUITMENT;

    uint32_t ss_id = 2;
    parameter ln_R0;
    parameter h;
    std::string name = "ricker";

    ricker() : recruitment_base() {
        recruitment_base::recruitment_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }
};

class beverton_holt : public recruitment_base {
public:
    uint32_t category = RECRUITMENT;
    uint32_t ss_id = 3;
    parameter ln_R0;
    parameter h;
    std::string name = "beverton_holt";

    beverton_holt() : recruitment_base() {
        recruitment_base::recruitment_objects[this->id] = this;
    }
    
    beverton_holt(const beverton_holt& other) :
    category(other.category), ss_id(other.ss_id), 
    ln_R0(other.ln_R0), h(other.h), name(other.name) {
    }


    uint32_t get_id() {
        return this->id;
    }
};

class hockey_stick : public recruitment_base {
public:
    uint32_t category = RECRUITMENT;
    uint32_t ss_id = 5;
    parameter ln_R0;
    parameter h;
    parameter Rmin;
    std::string name = "hockey_stick";

    hockey_stick() : recruitment_base() {
        recruitment_base::recruitment_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }
};

class selectivity_base : public model_base {
public:
    uint32_t category = SELECTIVITY;
    static std::map<uint32_t, selectivity_base*> selectivity_objects;
    static uint32_t id_g;
    uint32_t id;

    selectivity_base() {
        this->id = selectivity_base::id_g++;
        model_base::model_objects.push_back(this);
    }
};
std::map<uint32_t, selectivity_base*> selectivity_base::selectivity_objects;
uint32_t selectivity_base::id_g = 1;

class constant_selectivity : public selectivity_base {
public:
    uint32_t category = SELECTIVITY;
    uint32_t ss_id = 0;
    std::string name = "constant_selectivity";

    constant_selectivity() : selectivity_base() {
        selectivity_base::selectivity_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class logistic_selectivity : public selectivity_base {
public:
    uint32_t category = SELECTIVITY;
    uint32_t ss_id = 1;
    parameter slope;
    parameter median;
    std::string name = "logistic_selectivity";

    logistic_selectivity() : selectivity_base() {
        selectivity_base::selectivity_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class double_logistic_selectivity : public selectivity_base {
public:
    uint32_t category = SELECTIVITY;
    uint32_t ss_id = 8;
    parameter asc_slope;
    parameter asc_median;
    parameter desc_slope;
    parameter desc_median;
    std::string name = "double_logistic_selectivity";

    double_logistic_selectivity() : selectivity_base() {
        selectivity_base::selectivity_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class double_normal_selectivity : public selectivity_base {
public:
    uint32_t category = SELECTIVITY;
    uint32_t ss_id = 22;
    std::string name = "double_normal_selectivity";

    double_normal_selectivity() : selectivity_base() {
        selectivity_base::selectivity_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class growth_base : public model_base {
public:
    uint32_t category = GROWTH;
    static std::map<uint32_t, growth_base*> growth_objects;
    static uint32_t id_g;
    uint32_t id;

    growth_base() {
        this->id = growth_base::id_g++;
        model_base::model_objects.push_back(this);
    }
};
std::map<uint32_t, growth_base*> growth_base::growth_objects;
uint32_t growth_base::id_g = 1;

class von_bertalanffy : public growth_base {
public:
    uint32_t category = GROWTH;
    uint32_t ss_id = 1;
    std::string name = "von_bertalanffy";

    von_bertalanffy() : growth_base() {
        growth_base::growth_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class schnute : public growth_base {
public:
    uint32_t category = GROWTH;
    uint32_t ss_id = 2;
    std::string name = "schnute";

    schnute() : growth_base() {
        growth_base::growth_objects[this->id] = this;
    }

    uint32_t get_id() {
        return this->id;
    }

};

class ht4sa_ensemble {
public:

    Rcpp::IntegerVector recruitment_units;
    Rcpp::IntegerVector selectivity_units;
    Rcpp::IntegerVector growth_units;

    void add_recruitment_unit(uint32_t id) {
        recruitment_units.push_back(id);
    }

    void add_selectivity_unit(uint32_t id) {
        selectivity_units.push_back(id);
    }

    void add_growth_unit(uint32_t id) {
        growth_units.push_back(id);
    }

    Rcpp::IntegerVector get_recruitment_units() {
        return this->recruitment_units;
    }

    Rcpp::IntegerVector get_selectivity_units() {
        return this->selectivity_units;
    }

    Rcpp::IntegerVector get_growth_units() {
        return this->growth_units;
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
    Rcpp::StringVector fleetnames;
    std::string sourcefile;
    std::string type;
    std::string ReadVersion;
    bool eof;
    double EmpiricalWAA;
    double N_GP;
    double N_platoon;
    double recr_dist_method;
    double recr_global_area;
    double recr_dist_read;
    double recr_dist_inx;
    Rcpp::DataFrame recr_dist_pattern;
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

    Rcpp::class_<ricker>("ricker")
            .constructor()
            .field("category", &ricker::category)
            .field("ss_id", &ricker::ss_id)
            .field("name", &ricker::name)
            .field("ln_R0", &ricker::ln_R0)
            .field("h", &ricker::h)
            .method("get_id", &ricker::get_id);

    Rcpp::class_<beverton_holt>("beverton_holt")
            .constructor()
            .field("category", &beverton_holt::category)
            .field("ss_id", &beverton_holt::ss_id)
            .field("name", &beverton_holt::name)
            .field("ln_R0", &beverton_holt::ln_R0)
            .field("h", &beverton_holt::h)
            .method("get_id", &beverton_holt::get_id);

    Rcpp::class_<hockey_stick>("hockey_stick")
            .constructor()
            .field("category", &hockey_stick::category)
            .field("ss_id", &hockey_stick::ss_id)
            .field("name", &hockey_stick::name)
            .field("ln_R0", &hockey_stick::ln_R0)
            .field("h", &hockey_stick::h)
            .field("Rmin", &hockey_stick::Rmin)
            .method("get_id", &hockey_stick::get_id);

    Rcpp::class_<constant_selectivity>("constant_selectivity")
            .constructor()
            .field("category", &constant_selectivity::category)
            .field("ss_id", &constant_selectivity::ss_id)
            .field("name", &constant_selectivity::name)
            .method("get_id", &constant_selectivity::get_id);

    Rcpp::class_<logistic_selectivity>("logistic_selectivity")
            .constructor()
            .field("category", &logistic_selectivity::category)
            .field("ss_id", &logistic_selectivity::ss_id)
            .field("name", &logistic_selectivity::name)
            .field("median", &logistic_selectivity::median)
            .field("slope", &logistic_selectivity::slope)
            .method("get_id", &logistic_selectivity::get_id);

    Rcpp::class_<double_logistic_selectivity>("double_logistic_selectivity")
            .constructor()
            .field("category", &double_logistic_selectivity::category)
            .field("ss_id", &double_logistic_selectivity::ss_id)
            .field("name", &double_logistic_selectivity::name)
            .field("asc_median", &double_logistic_selectivity::asc_median)
            .field("asc_slope", &double_logistic_selectivity::asc_slope)
            .field("desc_median", &double_logistic_selectivity::desc_median)
            .field("desc_slope", &double_logistic_selectivity::desc_slope)
            .method("get_id", &double_logistic_selectivity::get_id);

    Rcpp::class_<double_normal_selectivity>("double_normal_selectivity")
            .constructor()
            .field("category", &double_normal_selectivity::category)
            .field("ss_id", &double_normal_selectivity::ss_id)
            .field("name", &double_normal_selectivity::name)
            .method("get_id", &double_normal_selectivity::get_id);

    Rcpp::class_<von_bertalanffy>("von_bertalanffy")
            .constructor()
            .field("category", &von_bertalanffy::category)
            .field("ss_id", &von_bertalanffy::ss_id)
            .field("name", &von_bertalanffy::name)
            .method("get_id", &von_bertalanffy::get_id);

    Rcpp::class_<schnute>("schnute")
            .constructor()
            .field("category", &schnute::category)
            .field("ss_id", &schnute::ss_id)
            .field("name", &schnute::name)
            .method("get_id", &schnute::get_id);

    Rcpp::class_<ht4sa_ensemble>("ht4sa_ensemble")
            .constructor()
            .method("add_recruitment_unit", &ht4sa_ensemble::add_recruitment_unit)
            .method("add_selectivity_unit", &ht4sa_ensemble::add_selectivity_unit)
            .method("add_growth_unit", &ht4sa_ensemble::add_growth_unit)
            .method("get_recruitment_units", &ht4sa_ensemble::get_recruitment_units)
            .method("get_selectivity_units", &ht4sa_ensemble::get_selectivity_units)
            .method("get_growth_units", &ht4sa_ensemble::get_growth_units);


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
            .field("fleetnames", &ht4sa_ss_control::fleetnames)
            .field("sourcefile", &ht4sa_ss_control::sourcefile)
            .field("type", &ht4sa_ss_control::type)
            .field("ReadVersion", &ht4sa_ss_control::ReadVersion)
            .field("eof", &ht4sa_ss_control::eof)
            .field("EmpiricalWAA", &ht4sa_ss_control::EmpiricalWAA)
            .field("N_GP", &ht4sa_ss_control::N_GP)
            .field("N_platoon", &ht4sa_ss_control::N_platoon)
            .field("recr_dist_method", &ht4sa_ss_control::recr_dist_method)
            .field("recr_global_area", &ht4sa_ss_control::recr_global_area)
            .field("recr_dist_read", &ht4sa_ss_control::recr_dist_read)
            .field("recr_dist_inx", &ht4sa_ss_control::recr_dist_inx)
            .field("recr_dist_pattern", &ht4sa_ss_control::recr_dist_pattern)
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
