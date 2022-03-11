
/*
 * File:   ht4sa.hpp
 * Authors:  Nicholas Ducharme-Barth, Megumi Oshima, Matthew Supernaw,
 *
 * Created on January 23, 2020, 9:03 AM
 */

#ifndef HT4SA_HPP
#define HT4SA_HPP

#include <cstdlib>
#include <iostream>

#include "mas/MAS.hpp"
#include "mas/Options.hpp"
#include "mas/ObjectiveFunction.hpp"

#include <RcppCommon.h>
#include <Rcpp.h>

#define  Rcout std::cout

using namespace Rcpp;

/**
 * Using the Google R style guide.
 */

template<typename T1, typename T2, typename T3>
struct triple {
    T1 first;
    T2 second;
    T3 third;

    triple(T1 first, T2 second, T3 third) :
    first(first), second(second), third(third) {
    }

};

struct Parameter {
public:
    double value;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();
    double lambda = 1.0;
    int phase = 1;
    bool estimated = false;

    Parameter(const Parameter &other) :
    value(other.value), min(other.min), max(other.max), phase(
    other.phase), estimated(other.estimated) {
    }

    Parameter(double value = 0.0) :
    value(value) {
    }

    virtual ~Parameter() {
    }

};

enum MASObjectType {
    SUBMODEL = 0, DATAOBJECT
};

class MASSubModel {
public:

    MASObjectType sub_model_type = SUBMODEL;

    typedef typename mas::VariableTrait<double>::variable variable;

    static std::vector<MASSubModel*> submodels;

    int references = 0;

    virtual ~MASSubModel() {
    }

    virtual void AddToMAS(mas::Information<double> &info) {

    }

    virtual void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

    }

    virtual void DataToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

    }

    void InitializeParameter(mas::ModelObject<double> *model, variable &v,
            const Parameter &p, const std::string &name) {
        v = p.value;
        if (p.estimated) {
            mas::VariableTrait<double>::SetName(v, name);
            if (p.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(v, p.min);
            }

            if (p.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(v, p.max);
            }
            model->Register(v, p.phase);
        }

    }

    static void GenerateArrayObject(rapidjson::Document &document,
            rapidjson::Value &array, const Rcpp::NumericVector &darray,
            int dimensions, size_t imax, size_t jmax, size_t kmax) {

        size_t i, j, k;

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();

        if (dimensions == 1) {

            for (i = 0; i < imax; i++) {
                array.PushBack(darray[i], allocator);
            }

        } else if (dimensions == 2) {
            for (i = 0; i < imax; i++) {
                array.PushBack(rapidjson::Value(rapidjson::kArrayType),
                        allocator);
                for (j = 0; j < jmax; j++) {
                    size_t index = i * jmax + j;
                    array[i].PushBack(darray[index], allocator);
                }
            }

        } else if (dimensions == 3) {

            for (i = 0; i < imax; i++) {
                array.PushBack(rapidjson::Value(rapidjson::kArrayType),
                        allocator);
                for (j = 0; j < jmax; j++) {
                    array[i].PushBack(rapidjson::Value(rapidjson::kArrayType),
                            allocator);
                    for (k = 0; k < kmax; k++) {

                        size_t index = i * jmax * kmax + j * kmax + k;
                        array[i][j].PushBack(darray[index], allocator);
                    }
                }
            }

        }

    }
};

std::vector<MASSubModel*> MASSubModel::submodels;

enum SexType {
    MALE = 0, FEMALE, UNDIFFERENTIATED, UNKOWN
};

SexType GetSexType(std::string sex) {
    if (sex == "males") {
        return MALE;
    } else if (sex == "females") {
        return FEMALE;
    } else if (sex == "undifferentiated") {
        return UNDIFFERENTIATED;
    } else {

        std::cout << "MAS Error: unknown sex type \"" << sex << "\".\n";

    }
    return UNDIFFERENTIATED;
}

/**
 * Selectivity
 */

class SelectivityBase : public MASSubModel {
protected:

public:
    Parameter sigma;
    Parameter sigma2;
    Parameter cv;
    static int id_g;

    SelectivityBase() {

    }

    virtual ~SelectivityBase() {
    }
};

int SelectivityBase::id_g = 1;

class LogisticSelectivity : public SelectivityBase {
public:

    Parameter a50; //age of 50% selectivity
    Parameter slope; //the rate of increase in selectivity at a50
    int id;
    Parameter sigma;
    Parameter sigma2;
    Parameter cv;

    LogisticSelectivity() : SelectivityBase() {
        cv.value = 0.05;
        this->id = SelectivityBase::id_g++;
        LogisticSelectivity::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    LogisticSelectivity(const LogisticSelectivity &other) :
    a50(other.a50), slope(other.slope), id(other.id) {
    }

    virtual ~LogisticSelectivity() {
    }

    double Evaluate(double a) {
        mas::VariableTrait<double>::variable A;
        mas::VariableTrait<double>::variable ret;
        mas::VariableTrait<double>::SetValue(A, a);
        atl::intrusive_ptr<mas::LogisticSel<double> > sel =
                new mas::LogisticSel<double>();
        mas::LogisticSel<double> *selex = sel.get();
        selex->a50 = this->a50.value;
        selex->s = this->slope.value;
        ret = selex->Evaluate(A);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::LogisticSel<double> > sel =
                new mas::LogisticSel<double>();


        mas::VariableTrait<double>::SetValue(sel->a50, a50.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->a50, a50.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->a50, a50.max);
        sel->id = id;

        mas::VariableTrait<double>::SetValue(sel->s, slope.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->s, slope.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->s, slope.max);


        mas::VariableTrait<double>::SetValue(sel->cv, this->cv.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->cv, this->cv.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->cv, this->cv.max);


        mas::VariableTrait<double>::SetValue(sel->sigma, sigma.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma, sigma.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma, sigma.max);


        mas::VariableTrait<double>::SetValue(sel->sigma2, sigma2.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma2, sigma2.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma2, sigma2.max);

        if (a50.estimated) {
            std::cout << "\n\nregistering a50\n";
            std::stringstream ss;
            ss << "logistic_selectivity_a50_" << id;
            mas::VariableTrait<double>::SetName(sel->a50, ss.str());
            sel->Register(sel->a50, a50.phase);
            sel->lambdas.push_back(a50.lambda);
        }

        sel->id = id;
        if (this->slope.estimated) {

            std::cout << "\n\nregistering slope\n";
            std::stringstream ss;
            ss << "logistic_selectivity_slope_" << id;
            mas::VariableTrait<double>::SetName(sel->s, ss.str());
            sel->Register(sel->s, slope.phase);
            sel->lambdas.push_back(slope.lambda);
        }
        info.selectivity_models[sel->id] = sel;

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "logistic", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value a50(rapidjson::kObjectType);
        a50.AddMember("value", this->a50.value, allocator);
        if (this->a50.estimated) {
            a50.AddMember("estimated", "true", allocator);
            a50.AddMember("min", this->a50.min, allocator);
            a50.AddMember("max", this->a50.max, allocator);
            a50.AddMember("phase", this->a50.phase, allocator);
        } else {
            a50.AddMember("estimated", "false", allocator);
            a50.AddMember("min", this->a50.min, allocator);
            a50.AddMember("max", this->a50.max, allocator);
            a50.AddMember("phase", this->a50.phase, allocator);
        }

        parameters.AddMember("a50", a50, allocator);

        rapidjson::Value s(rapidjson::kObjectType);
        s.AddMember("value", this->slope.value, allocator);
        if (this->slope.estimated) {

            s.AddMember("estimated", "true", allocator);
            s.AddMember("min", this->slope.min, allocator);
            s.AddMember("max", this->slope.max, allocator);
            s.AddMember("phase", this->slope.phase, allocator);
        } else {
            s.AddMember("estimated", "false", allocator);
            s.AddMember("min", this->slope.min, allocator);
            s.AddMember("max", this->slope.max, allocator);
            s.AddMember("phase", this->slope.phase, allocator);
        }

        parameters.AddMember("s", s, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        document.AddMember("selectivity", selectivity, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "logistic", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value a50(rapidjson::kObjectType);
        a50.AddMember("value", this->a50.value, allocator);
        if (this->a50.estimated) {
            a50.AddMember("estimated", "true", allocator);
            a50.AddMember("min", this->a50.min, allocator);
            a50.AddMember("max", this->a50.max, allocator);
            a50.AddMember("phase", this->a50.phase, allocator);
        } else {
            a50.AddMember("estimated", "false", allocator);
            a50.AddMember("min", this->a50.min, allocator);
            a50.AddMember("max", this->a50.max, allocator);
            a50.AddMember("phase", this->a50.phase, allocator);
        }

        parameters.AddMember("a50", a50, allocator);

        rapidjson::Value s(rapidjson::kObjectType);
        s.AddMember("value", this->slope.value, allocator);
        if (this->slope.estimated) {

            s.AddMember("estimated", "true", allocator);
            s.AddMember("min", this->slope.min, allocator);
            s.AddMember("max", this->slope.max, allocator);
            s.AddMember("phase", this->slope.phase, allocator);
        } else {
            s.AddMember("estimated", "false", allocator);
            s.AddMember("min", this->slope.min, allocator);
            s.AddMember("max", this->slope.max, allocator);
            s.AddMember("phase", this->slope.phase, allocator);
        }

        parameters.AddMember("s", s, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        selex.PushBack(selectivity, allocator);
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::selectivity_model_iterator it;
        it = info.selectivity_models.find(this->id);
        if (it != info.selectivity_models.end()) {

            atl::intrusive_ptr<mas::SelectivityBase<double> > sel = (*it).second;
            mas::LogisticSel<double> *lsel =
                    (mas::LogisticSel<double>*) sel.get();
            this->a50.value = lsel->a50.GetValue();
            this->slope.value = lsel->s.GetValue();
        }
    }

    static std::map<int, LogisticSelectivity*> initialized_models;
    typedef typename std::map<int, LogisticSelectivity*>::iterator model_iterator;
};

std::map<int, LogisticSelectivity*> LogisticSelectivity::initialized_models;

class DoubleLogisticSelectivity : public SelectivityBase {
public:

    int id;
    Parameter alpha_asc; //ascending alpha
    Parameter beta_asc; // ascending beta
    Parameter alpha_desc; // descending alpha
    Parameter beta_desc; // descending beta
    Parameter sigma;
    Parameter sigma2;
    Parameter cv;

    DoubleLogisticSelectivity() : SelectivityBase() {
        cv.value = 0.05;

        this->id = SelectivityBase::id_g++;
        DoubleLogisticSelectivity::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);

    }

    DoubleLogisticSelectivity(const DoubleLogisticSelectivity &other) :
    id(other.id), alpha_asc(other.alpha_asc), beta_asc(other.beta_asc), alpha_desc(
    other.alpha_desc), beta_desc(other.beta_desc) {
    }

    virtual ~DoubleLogisticSelectivity() {
    }

    double Evaluate(double a) {
        mas::VariableTrait<double>::variable A;
        mas::VariableTrait<double>::variable ret;
        mas::VariableTrait<double>::SetValue(A, a);
        atl::intrusive_ptr<mas::DoubleLogisticSel<double> > sel =
                new mas::DoubleLogisticSel<double>();
        mas::DoubleLogisticSel<double> *selex = sel.get();
        selex->alpha_asc = this->alpha_asc.value;
        selex->alpha_desc = this->alpha_desc.value;
        selex->beta_asc = this->beta_asc.value;
        selex->beta_desc = this->beta_desc.value;
        ret = selex->Evaluate(A);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::DoubleLogisticSel<double> > sel =
                new mas::DoubleLogisticSel<double>();
        sel->id = id;
        mas::DoubleLogisticSel<double> *selex = sel.get();

        mas::VariableTrait<double>::SetValue(sel->cv, this->cv.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->cv, this->cv.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->cv, this->cv.max);

        mas::VariableTrait<double>::SetValue(sel->sigma, sigma.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma, sigma.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma, sigma.max);


        mas::VariableTrait<double>::SetValue(sel->sigma2, sigma2.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma2, sigma2.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma2, sigma2.max);

        mas::VariableTrait<double>::SetValue(selex->alpha_asc,
                this->alpha_asc.value);
        mas::VariableTrait<double>::SetMinBoundary(selex->alpha_asc,
                this->alpha_asc.min);
        mas::VariableTrait<double>::SetMaxBoundary(selex->alpha_asc,
                this->alpha_asc.max);
        if (this->alpha_asc.estimated) {
            std::stringstream ss;
            ss << "double_logistic_selectivity_alpha_asc_" << id;
            mas::VariableTrait<double>::SetName(selex->alpha_asc, ss.str());
            sel->Register(selex->alpha_asc, this->alpha_asc.phase);
            sel->lambdas.push_back(alpha_asc.lambda);
        }

        mas::VariableTrait<double>::SetValue(selex->beta_asc,
                this->beta_asc.value);
        mas::VariableTrait<double>::SetMinBoundary(selex->beta_asc,
                this->beta_asc.min);
        mas::VariableTrait<double>::SetMaxBoundary(selex->beta_asc,
                this->beta_asc.max);
        if (this->alpha_asc.estimated) {
            std::stringstream ss;
            ss << "double_logistic_selectivity_beta_asc_" << id;
            mas::VariableTrait<double>::SetName(selex->beta_asc, ss.str());
            sel->Register(selex->beta_asc, this->beta_asc.phase);
            sel->lambdas.push_back(beta_asc.lambda);
        }

        mas::VariableTrait<double>::SetValue(selex->alpha_desc,
                this->alpha_desc.value);
        mas::VariableTrait<double>::SetMinBoundary(selex->alpha_desc,
                this->alpha_desc.min);
        mas::VariableTrait<double>::SetMaxBoundary(selex->alpha_desc,
                this->alpha_desc.max);
        if (this->alpha_asc.estimated) {
            std::stringstream ss;
            ss << "double_logistic_selectivity_alpha_desc_" << id;
            mas::VariableTrait<double>::SetName(selex->alpha_desc, ss.str());
            sel->Register(selex->alpha_desc, this->alpha_desc.phase);
            sel->lambdas.push_back(alpha_desc.lambda);
        }

        mas::VariableTrait<double>::SetValue(selex->beta_desc,
                this->beta_desc.value);
        mas::VariableTrait<double>::SetMinBoundary(selex->beta_desc,
                this->beta_desc.min);
        mas::VariableTrait<double>::SetMaxBoundary(selex->beta_desc,
                this->beta_desc.max);
        if (this->alpha_asc.estimated) {

            std::stringstream ss;
            ss << "double_logistic_selectivity_beta_desc_" << id;
            mas::VariableTrait<double>::SetName(selex->beta_desc, ss.str());
            sel->Register(selex->beta_desc, this->beta_desc.phase);
            sel->lambdas.push_back(beta_desc.lambda);
        }

        info.selectivity_models[sel->id] = sel;

    }

    void ExtractFromMAS(mas::Information<double> &info) {

        typename mas::Information<double>::selectivity_model_iterator it;
        it = info.selectivity_models.find(this->id);
        if (it != info.selectivity_models.end()) {

            atl::intrusive_ptr<mas::SelectivityBase<double> > sel = (*it).second;
            mas::DoubleLogisticSel<double> *lsel = (mas::DoubleLogisticSel<
                    double>*) sel.get();
            this->alpha_asc.value = lsel->alpha_asc.GetValue();
            this->beta_asc.value = lsel->beta_asc.GetValue();
            this->alpha_desc = lsel->alpha_desc.GetValue();
            this->beta_desc = lsel->beta_desc.GetValue();
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "double_logistic", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value alpha_asc(rapidjson::kObjectType);
        alpha_asc.AddMember("value", this->alpha_asc.value, allocator);
        if (this->alpha_asc.estimated) {
            alpha_asc.AddMember("estimated", "true", allocator);
            alpha_asc.AddMember("min", this->alpha_asc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_asc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_asc.phase, allocator);
        } else {
            alpha_asc.AddMember("estimated", "false", allocator);
            alpha_asc.AddMember("min", this->alpha_asc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_asc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_asc.phase, allocator);
        }

        parameters.AddMember("alpha_asc", alpha_asc, allocator);

        rapidjson::Value beta_asc(rapidjson::kObjectType);
        beta_asc.AddMember("value", this->beta_asc.value, allocator);
        if (this->beta_asc.estimated) {

            beta_asc.AddMember("estimated", "true", allocator);
            beta_asc.AddMember("min", this->beta_asc.min, allocator);
            beta_asc.AddMember("max", this->beta_asc.max, allocator);
            beta_asc.AddMember("phase", this->beta_asc.phase, allocator);
        } else {

            beta_asc.AddMember("estimated", "false", allocator);
            beta_asc.AddMember("min", this->beta_asc.min, allocator);
            beta_asc.AddMember("max", this->beta_asc.max, allocator);
            beta_asc.AddMember("phase", this->beta_asc.phase, allocator);
        }

        parameters.AddMember("beta_asc", beta_asc, allocator);

        rapidjson::Value alpha_desc(rapidjson::kObjectType);
        alpha_desc.AddMember("value", this->alpha_desc.value, allocator);
        if (this->alpha_desc.estimated) {
            alpha_asc.AddMember("estimated", "true", allocator);
            alpha_asc.AddMember("min", this->alpha_desc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_desc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_desc.phase, allocator);
        } else {
            alpha_asc.AddMember("estimated", "false", allocator);
            alpha_asc.AddMember("min", this->alpha_desc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_desc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_desc.phase, allocator);
        }

        parameters.AddMember("alpha_desc", alpha_desc, allocator);

        rapidjson::Value beta_desc(rapidjson::kObjectType);
        beta_desc.AddMember("value", this->beta_desc.value, allocator);
        if (this->beta_desc.estimated) {

            beta_desc.AddMember("estimated", "true", allocator);
            beta_desc.AddMember("min", this->beta_desc.min, allocator);
            beta_desc.AddMember("max", this->beta_desc.max, allocator);
            beta_desc.AddMember("phase", this->beta_desc.phase, allocator);
        } else {

            beta_desc.AddMember("estimated", "false", allocator);
            beta_desc.AddMember("min", this->beta_desc.min, allocator);
            beta_desc.AddMember("max", this->beta_desc.max, allocator);
            beta_desc.AddMember("phase", this->beta_desc.phase, allocator);
        }

        parameters.AddMember("beta_desc", beta_desc, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        document.AddMember("selectivity", selectivity, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "double_logistic", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value alpha_asc(rapidjson::kObjectType);
        alpha_asc.AddMember("value", this->alpha_asc.value, allocator);
        if (this->alpha_asc.estimated) {
            alpha_asc.AddMember("estimated", "true", allocator);
            alpha_asc.AddMember("min", this->alpha_asc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_asc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_asc.phase, allocator);
        } else {
            alpha_asc.AddMember("estimated", "false", allocator);
            alpha_asc.AddMember("min", this->alpha_asc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_asc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_asc.phase, allocator);
        }

        parameters.AddMember("alpha_asc", alpha_asc, allocator);

        rapidjson::Value beta_asc(rapidjson::kObjectType);
        beta_asc.AddMember("value", this->beta_asc.value, allocator);
        if (this->beta_asc.estimated) {

            beta_asc.AddMember("estimated", "true", allocator);
            beta_asc.AddMember("min", this->beta_asc.min, allocator);
            beta_asc.AddMember("max", this->beta_asc.max, allocator);
            beta_asc.AddMember("phase", this->beta_asc.phase, allocator);
        } else {

            beta_asc.AddMember("estimated", "false", allocator);
            beta_asc.AddMember("min", this->beta_asc.min, allocator);
            beta_asc.AddMember("max", this->beta_asc.max, allocator);
            beta_asc.AddMember("phase", this->beta_asc.phase, allocator);
        }

        parameters.AddMember("beta_asc", beta_asc, allocator);

        rapidjson::Value alpha_desc(rapidjson::kObjectType);
        alpha_desc.AddMember("value", this->alpha_desc.value, allocator);
        if (this->alpha_desc.estimated) {
            alpha_asc.AddMember("estimated", "true", allocator);
            alpha_asc.AddMember("min", this->alpha_desc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_desc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_desc.phase, allocator);
        } else {
            alpha_asc.AddMember("estimated", "false", allocator);
            alpha_asc.AddMember("min", this->alpha_desc.min, allocator);
            alpha_asc.AddMember("max", this->alpha_desc.max, allocator);
            alpha_asc.AddMember("phase", this->alpha_desc.phase, allocator);
        }

        parameters.AddMember("alpha_desc", alpha_desc, allocator);

        rapidjson::Value beta_desc(rapidjson::kObjectType);
        beta_desc.AddMember("value", this->beta_desc.value, allocator);
        if (this->beta_desc.estimated) {

            beta_desc.AddMember("estimated", "true", allocator);
            beta_desc.AddMember("min", this->beta_desc.min, allocator);
            beta_desc.AddMember("max", this->beta_desc.max, allocator);
            beta_desc.AddMember("phase", this->beta_desc.phase, allocator);
        } else {

            beta_desc.AddMember("estimated", "false", allocator);
            beta_desc.AddMember("min", this->beta_desc.min, allocator);
            beta_desc.AddMember("max", this->beta_desc.max, allocator);
            beta_desc.AddMember("phase", this->beta_desc.phase, allocator);
        }

        parameters.AddMember("beta_desc", beta_desc, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        selex.PushBack(selectivity, allocator);
    }

    static std::map<int, DoubleLogisticSelectivity*> initialized_models;
    typedef typename std::map<int, DoubleLogisticSelectivity*>::iterator model_iterator;
};

std::map<int, DoubleLogisticSelectivity*> DoubleLogisticSelectivity::initialized_models;

class AgeBasedSelectivity : public SelectivityBase {
public:

    Rcpp::NumericVector values;
    Rcpp::IntegerVector estimate_age;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();
    bool estimated = false;
    int phase = 1;
    int id;
    Parameter sigma;
    Parameter sigma2;
    Parameter cv;
    
    AgeBasedSelectivity() : SelectivityBase() {
        cv.value = 0.05;
        this->id = SelectivityBase::id_g++;
        AgeBasedSelectivity::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~AgeBasedSelectivity() {
    }

    double Evaluate(double a) {
        mas::VariableTrait<double>::variable A;
        mas::VariableTrait<double>::variable ret;
        mas::VariableTrait<double>::SetValue(A, a);
        atl::intrusive_ptr<mas::AgeBased<double> > sel = new mas::AgeBased<
                double>();
        mas::AgeBased<double> *selex = sel.get();
        for (int i = 0; i < this->values.size(); i++) {
            selex->w.push_back(this->values[i]);
        }

        ret = selex->Evaluate(A);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        if (this->values.size() != info.ages.size()) {
            info.valid_configuration = false;
            std::cout << "mismatch in nages and AgeBasedSelectivity vecrtor\n";
            return;
        }
        atl::intrusive_ptr<mas::AgeBased<double> > sel = new mas::AgeBased<
                double>();

        mas::VariableTrait<double>::SetValue(sel->cv, this->cv.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->cv, this->cv.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->cv, this->cv.max);

        mas::VariableTrait<double>::SetValue(sel->sigma, sigma.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma, sigma.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma, sigma.max);


        mas::VariableTrait<double>::SetValue(sel->sigma2, sigma2.value);
        mas::VariableTrait<double>::SetMinBoundary(sel->sigma2, sigma2.min);
        mas::VariableTrait<double>::SetMaxBoundary(sel->sigma2, sigma2.max);

        mas::AgeBased<double> *selex = sel.get();

        selex->id = this->id;
        for (int i = 0; i < this->values.size(); i++) {
            selex->w.push_back(this->values[i]);
        }

        if (this->estimated) {
            if (this->estimate_age.size() != this->values.size()
                    || this->estimate_age.size() == 0) {

                if (this->estimate_age.size() > 0) {
                    std::cout
                            << "Warning: Vector \"estimate_age\" for age based "
                            "selectivity model \"" << this->id
                            << "\" not 0 or "
                            "values.size(). Resizing and setting all values to 1.\n";
                    mas::mas_log
                            << "Warning: Vector \"estimate_age\" for age based "
                            "selectivity model \"" << this->id
                            << "\" not 0 or "
                            "values.size(). Resizing and setting all values to 1.\n";
                }

                for (int i = 0; i < this->values.size(); i++) {
                    this->estimate_age.push_back(1);
                }
            }

            for (int i = 0; i < this->values.size(); i++) {

                std::stringstream ss;
                ss << "age_base_selectivity[" << i << "]_" << this->id;
                selex->w[i].SetName(ss.str());
                selex->selex[info.ages[i]] = this->values[i];
                selex->selex[info.ages[i]].SetName(ss.str());
                if (this->estimate_age[i] > 0) {
                    selex->Register(selex->selex[info.ages[i]], this->phase);
                }
            }
        }
        info.selectivity_models[selex->id] = sel;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::selectivity_model_iterator it;
        it = info.selectivity_models.find(this->id);
        if (it != info.selectivity_models.end()) {
            atl::intrusive_ptr<mas::SelectivityBase<double> > sel = (*it).second;
            mas::AgeBased<double> *lsel = (mas::AgeBased<double>*) sel.get();
            for (int i = 0; i < lsel->w.size(); i++) {

                this->values[i] = lsel->w[i].GetValue();
            }
        }

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "age_based", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value s(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        rapidjson::Value min(rapidjson::kArrayType);
        rapidjson::Value max(rapidjson::kArrayType);
        for (int i = 0; i < values.size(); i++) {
            vals.PushBack(this->values[i], allocator);
            min.PushBack(this->min, allocator);
            max.PushBack(this->max, allocator);
        }
        s.AddMember("values", vals, allocator);
        s.AddMember("min", min, allocator);
        s.AddMember("max", max, allocator);
        if (this->estimated) {
            s.AddMember("estimated", "true", allocator);
            s.AddMember("min", this->min, allocator);
            s.AddMember("max", this->max, allocator);
            s.AddMember("phase", this->phase, allocator);
        } else {
            s.AddMember("estimated", "false", allocator);
            s.AddMember("min", this->min, allocator);
            s.AddMember("max", this->max, allocator);
            s.AddMember("phase", this->phase, allocator);
        }

        parameters.AddMember("s", s, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        document.AddMember("selectivity", selectivity, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value selectivity(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        selectivity.AddMember("model", "age_based", allocator);
        selectivity.AddMember("id", this->id, allocator);
        rapidjson::Value s(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        rapidjson::Value min(rapidjson::kArrayType);
        rapidjson::Value max(rapidjson::kArrayType);
        for (int i = 0; i < values.size(); i++) {
            vals.PushBack(this->values[i], allocator);
            min.PushBack(this->min, allocator);
            max.PushBack(this->max, allocator);
        }
        s.AddMember("values", vals, allocator);
        s.AddMember("min", min, allocator);
        s.AddMember("max", max, allocator);
        if (this->estimated) {
            s.AddMember("estimated", "true", allocator);
            s.AddMember("min", this->min, allocator);
            s.AddMember("max", this->max, allocator);
            s.AddMember("phase", this->phase, allocator);
        } else {
            s.AddMember("estimated", "false", allocator);
            s.AddMember("min", this->min, allocator);
            s.AddMember("max", this->max, allocator);
            s.AddMember("phase", this->phase, allocator);
        }

        parameters.AddMember("s", s, allocator);
        selectivity.AddMember("parameters", parameters, allocator);
        selex.PushBack(selectivity, allocator);
    }

    static std::map<int, AgeBasedSelectivity*> initialized_models;
    typedef typename std::map<int, AgeBasedSelectivity*>::iterator model_iterator;
};

std::map<int, AgeBasedSelectivity*> AgeBasedSelectivity::initialized_models;

/**
 * Mortality
 */
class FishingMortality : public MASSubModel {
public:
    static int id_g;
    int id;
    Rcpp::NumericVector values;
    bool estimate = false;
    int phase = 1;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();

    FishingMortality() {

        this->id = FishingMortality::id_g++;
        FishingMortality::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);

    }

    FishingMortality(const FishingMortality &other) :
    id(other.id), values(other.values), estimate(other.estimate), phase(
    other.phase), min(other.min), max(other.max) {

        std::cout << "Fishing Mortaility copy constructor";
    }

    virtual ~FishingMortality() {
    }

    void SetValues(Rcpp::NumericVector values) {

        this->values = values;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::FishingMortality<double> > fm =
                new mas::FishingMortality<double>();
        mas::FishingMortality<double> *m = fm.get();
        m->fishing_mortality_type = mas::ESTIMATED;
        fm->id = this->id;
        if (this->values.size() != info.nyears * info.nseasons) {
            std::cout
                    << "MAS Error: FishingMortality vector not equal to (nyears*nseasons)...";
            std::cout << this->values.size() << "!="
                    << (info.nyears * info.nseasons) << "\n";

            info.valid_configuration = false;
            return;
        }

        fm->fishing_mortality.resize(info.nyears);
        typedef typename mas::VariableTrait<double>::variable variable;
        int i = 0;
        for (int y = 0; y < info.nyears; y++) {
            fm->fishing_mortality[y].resize(info.nseasons);
            for (int s = 0; s < info.nseasons; s++) {
                fm->fishing_mortality[y][s] = variable(this->values[i++]);
            }
        }

        if (this->estimate) {
            for (int y = 0; y < info.nyears; y++) {
                fm->fishing_mortality[y].resize(info.nseasons);
                for (int s = 0; s < info.nseasons; s++) {
                    std::stringstream ss;
                    ss << "fishing_mortality[" << y << "][" << s << "]_"
                            << this->id;
                    mas::VariableTrait<double>::SetName(
                            fm->fishing_mortality[y][s], ss.str());
                    if (this->min != std::numeric_limits<double>::min()) {
                        mas::VariableTrait<double>::SetMinBoundary(
                                fm->fishing_mortality[y][s], this->min);
                    }
                    if (this->max != std::numeric_limits<double>::max()) {

                        mas::VariableTrait<double>::SetMaxBoundary(
                                fm->fishing_mortality[y][s], this->max);
                    }
                    std::cout << "Registering " << ss.str() << " "
                            << this->phase << "\n";
                    fm->Register(fm->fishing_mortality[y][s], this->phase);
                }
            }
        }

        info.fishing_mortality_models[m->id] = fm;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::fishing_mortality_model_iterator fit;
        fit = info.fishing_mortality_models.find(this->id);
        if (fit != info.fishing_mortality_models.end()) {
            int i = 0;
            for (int y = 0; y < info.nyears; y++) {
                for (int s = 0; s < info.nseasons; s++) {

                    this->values[i++] =
                            (*fit).second->fishing_mortality[y][s].GetValue();
                }
            }
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value fishing_mort(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->values.size(); i++) {
            rapidjson::Value entry(rapidjson::kArrayType);
            entry.PushBack(this->values[i], allocator);
            vals.PushBack(entry, allocator);
        }
        fishing_mort.AddMember("id", this->id, allocator);
        if (this->estimate) {
            parameters.AddMember("estimated", "true", allocator);
        } else {
            parameters.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("phase", this->phase, allocator);
        parameters.AddMember("min", this->min, allocator);
        parameters.AddMember("max", this->max, allocator);
        parameters.AddMember("values", vals, allocator);
        fishing_mort.AddMember("parameters", parameters, allocator);
        document.AddMember("fishing_mortality", fishing_mort, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &fm, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value fishing_mort(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->values.size(); i++) {
            rapidjson::Value entry(rapidjson::kArrayType);
            entry.PushBack(this->values[i], allocator);
            vals.PushBack(entry, allocator);
        }
        fishing_mort.AddMember("id", this->id, allocator);
        if (this->estimate) {
            parameters.AddMember("estimated", "true", allocator);
        } else {
            parameters.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("phase", this->phase, allocator);
        parameters.AddMember("min", this->min, allocator);
        parameters.AddMember("max", this->max, allocator);
        parameters.AddMember("values", vals, allocator);
        fishing_mort.AddMember("parameters", parameters, allocator);
        fm.PushBack(fishing_mort, allocator);
    }

    static std::map<int, FishingMortality*> initialized_models;
    typedef typename std::map<int, FishingMortality*>::iterator model_iterator;
};

std::map<int, FishingMortality*> FishingMortality::initialized_models;
int FishingMortality::id_g = 1;

class NaturalMortality : public MASSubModel {
public:
    static int id_g;
    int id;
    Rcpp::NumericVector values;
    bool estimate = false;
    int phase = 1;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();

    NaturalMortality() {

        this->id = NaturalMortality::id_g++;
        NaturalMortality::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);

    }

    virtual ~NaturalMortality() {
    }

    void SetValues(Rcpp::NumericVector values) {

        this->values = values;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        typedef typename mas::VariableTrait<double>::variable variable;

        atl::intrusive_ptr<mas::NaturalMortality<double> > nm =
                new mas::NaturalMortality<double>();
        mas::NaturalMortality<double> *m = nm.get();
        m->mortality_vector.resize(info.ages.size());
        m->id = this->id;

        for (int i = 0; i < info.ages.size(); i++) {
            m->mortality_vector[i] = variable(this->values[i]);
        }

        if (this->estimate) {
            for (int i = 0; i < info.ages.size(); i++) {
                std::stringstream ss;
                ss << "natural_mortality[" << i << "]_" << this->id;
                mas::VariableTrait<double>::SetName(m->mortality_vector[i],
                        ss.str());
                if (min != std::numeric_limits<double>::min()) {
                    mas::VariableTrait<double>::SetMinBoundary(
                            m->mortality_vector[i], min);
                }
                if (max != std::numeric_limits<double>::max()) {

                    mas::VariableTrait<double>::SetMaxBoundary(
                            m->mortality_vector[i], max);
                }
                m->Register(m->mortality_vector[i], phase);
            }
        }

        info.natural_mortality_models[this->id] = nm;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::natural_mortality_model_iterator nit;
        nit = info.natural_mortality_models.find(this->id);
        if (nit != info.natural_mortality_models.end()) {
            atl::intrusive_ptr<mas::NaturalMortality<double> > nm =
                    (*nit).second;
            mas::NaturalMortality<double> *m = nm.get();
            for (int i = 0; i < info.ages.size(); i++) {

                this->values[i] = m->mortality_vector[i].GetValue();
            }
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value natural_mort(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        rapidjson::Value min(rapidjson::kArrayType);
        rapidjson::Value max(rapidjson::kArrayType);
        for (int i = 0; i < this->values.size(); i++) {
            vals.PushBack(this->values[i], allocator);
            min.PushBack(this->min, allocator);
            max.PushBack(this->max, allocator);
        }
        natural_mort.AddMember("id", this->id, allocator);
        if (this->estimate) {
            parameters.AddMember("estimated", "true", allocator);
        } else {
            parameters.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("phase", this->phase, allocator);
        parameters.AddMember("values", vals, allocator);
        parameters.AddMember("min", min, allocator);
        parameters.AddMember("max", max, allocator);

        natural_mort.AddMember("parameters", parameters, allocator);
        document.AddMember("natural_mortality", natural_mort, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nm, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value natural_mort(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        rapidjson::Value min(rapidjson::kArrayType);
        rapidjson::Value max(rapidjson::kArrayType);
        for (int i = 0; i < this->values.size(); i++) {
            vals.PushBack(this->values[i], allocator);
            min.PushBack(this->min, allocator);
            max.PushBack(this->max, allocator);
        }
        natural_mort.AddMember("id", this->id, allocator);
        if (this->estimate) {
            parameters.AddMember("estimated", "true", allocator);
        } else {
            parameters.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("phase", this->phase, allocator);
        parameters.AddMember("values", vals, allocator);
        parameters.AddMember("min", min, allocator);
        parameters.AddMember("max", max, allocator);

        natural_mort.AddMember("parameters", parameters, allocator);
        nm.PushBack(natural_mort, allocator);
    }

    static std::map<int, NaturalMortality*> initialized_models;
    typedef typename std::map<int, NaturalMortality*>::iterator model_iterator;
};

std::map<int, NaturalMortality*> NaturalMortality::initialized_models;
int NaturalMortality::id_g = 1;

class InitialDeviations {
public:
    static int id_g;
    int id;
    Rcpp::NumericVector values;
    bool estimate = false;
    int phase = 1;
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();

    InitialDeviations() {
        this->id = InitialDeviations::id_g++;
        InitialDeviations::initialized_models[this->id] = this;
    }

    virtual ~InitialDeviations() {
    }

    void SetValues(Rcpp::NumericVector values) {

        this->values = values;
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &init_devs, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value initial_devs(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        rapidjson::Value min(rapidjson::kArrayType);
        rapidjson::Value max(rapidjson::kArrayType);
        for (int i = 0; i < this->values.size(); i++) {
            vals.PushBack(this->values[i], allocator);
            min.PushBack(this->min, allocator);
            max.PushBack(this->max, allocator);
        }
        initial_devs.AddMember("id", this->id, allocator);
        if (this->estimate) {
            parameters.AddMember("estimated", "true", allocator);
        } else {
            parameters.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("phase", this->phase, allocator);
        parameters.AddMember("values", vals, allocator);
        parameters.AddMember("min", min, allocator);
        parameters.AddMember("max", max, allocator);

        initial_devs.AddMember("parameters", parameters, allocator);
        init_devs.PushBack(initial_devs, allocator);

    }

    static std::map<int, InitialDeviations*> initialized_models;
    typedef typename std::map<int, InitialDeviations*>::iterator model_iterator;
};

std::map<int, InitialDeviations*> InitialDeviations::initialized_models;
int InitialDeviations::id_g = 1;

/**
 * Recruitment
 */
class RecruitmentBase : public MASSubModel {
protected:

public:
    static int id_g;

    virtual ~RecruitmentBase() {
    }
    Parameter sigma_r;

};
int RecruitmentBase::id_g = 1;

class RickerRecruitment : public RecruitmentBase {
public:

    Rcpp::NumericVector deviations;
    double deviations_min = std::numeric_limits<double>::min();
    double deviations_max = std::numeric_limits<double>::max();

    bool estimate_deviations = true;
    int deviation_phase = 1;
    bool constrained_deviations = true;

    Parameter R0;
    Parameter alpha;
    Parameter beta;
    bool use_bias_correction = false;
    int id;

    RickerRecruitment() {

        this->id = RecruitmentBase::id_g++;
        RickerRecruitment::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);

    }

    virtual ~RickerRecruitment() {
    }

    void SetDeviations(Rcpp::NumericVector values) {

        this->deviations = values;
    }

    double Evaluate(double SB0, double sb) {
        typedef typename mas::VariableTrait<double>::variable variable;
        variable ret;
        variable SB0_;
        variable sb_;

        mas::VariableTrait<double>::SetValue(SB0_, SB0);
        mas::VariableTrait<double>::SetValue(sb_, sb);

        atl::intrusive_ptr<mas::Ricker<double> > rec =
                new mas::Ricker<double>();

        mas::Ricker<double> *r = rec.get();
        r->SB0[1][1] = SB0_;
        r->R0 = this->R0.value;
        r->log_R0 = std::log(this->R0.value);
        r->alpha = this->alpha.value;
        r->beta = this->beta.value;
        ret = r->Evaluate(1, 1, sb_);

        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        typedef typename mas::VariableTrait<double>::variable variable;

        atl::intrusive_ptr<mas::Ricker<double> > rec =
                new mas::Ricker<double>();
        mas::Ricker<double> *r = rec.get();
        r->id = this->id;

        r->log_R0 = std::log(this->R0.value);


        r->use_bias_correction = this->use_bias_correction;

        if (this->R0.estimated) {
            std::stringstream ss;
            ss << "log_R0" << this->id;
            mas::VariableTrait<double>::SetName(r->log_R0, ss.str());
            if (this->R0.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->log_R0,
                        std::log(this->R0.min));
            }
            if (this->R0.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->log_R0,
                        std::log(this->R0.max));
            }
            r->Register(r->log_R0, this->R0.phase);
        }

        r->alpha = this->alpha.value;
        if (this->alpha.estimated) {
            std::stringstream ss;
            ss << "alpha" << this->id;
            mas::VariableTrait<double>::SetName(r->alpha, ss.str());
            if (this->alpha.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->alpha, alpha.min);
            }
            if (this->alpha.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->alpha, alpha.max);
            }
            r->Register(r->alpha, this->alpha.phase);
        }

        r->beta = this->beta.value;
        if (this->beta.estimated) {
            std::stringstream ss;
            ss << "beta" << this->id;
            mas::VariableTrait<double>::SetName(r->beta, ss.str());
            if (this->beta.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->beta, beta.min);
            }
            if (this->beta.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->beta, beta.max);
            }
            r->Register(r->beta, this->beta.phase);
        }

        r->recruitment_deviations.resize(this->deviations.size());
        for (int i = 0; i < this->deviations.size(); i++) {
            r->recruitment_deviations[i] = variable(this->deviations[i]);
        }
        r->recruitment_deviations_constrained = this->constrained_deviations;

        if (this->estimate_deviations) {
            r->estimating_recruitment_deviations = true;
            for (int i = 0; i < info.nyears; i++) {
                std::stringstream ss;
                ss << "recruitment_deviations[" << i << "]_" << this->id;
                mas::VariableTrait<double>::SetName(
                        r->recruitment_deviations[i], ss.str());
                if (this->deviations_min
                        != std::numeric_limits<double>::min()) {
                    mas::VariableTrait<double>::SetMinBoundary(
                            r->recruitment_deviations[i], this->deviations_min);
                }
                if (this->deviations_max
                        != std::numeric_limits<double>::max()) {

                    mas::VariableTrait<double>::SetMaxBoundary(
                            r->recruitment_deviations[i], this->deviations_max);
                }
                r->recruitment_deviations[i] = variable(this->deviations[i]);
                r->Register(r->recruitment_deviations[i],
                        this->deviation_phase);

            }
        }

        info.recruitment_models[this->id] = rec;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::recruitment_model_iterator rit;
        rit = info.recruitment_models.find(this->id);
        if (rit != info.recruitment_models.end()) {

            atl::intrusive_ptr<mas::RecruitmentBase<double> > rec =
                    (*rit).second;
            mas::Ricker<double> *r = (mas::Ricker<double>*) rec.get();
            this->R0.value = std::exp(r->log_R0.GetValue());
            this->alpha.value = r->alpha.GetValue();
            this->beta.value = r->beta.GetValue();
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "ricker", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value alpha(rapidjson::kObjectType);
        alpha.AddMember("value", this->alpha.value, allocator);
        if (this->alpha.estimated) {
            alpha.AddMember("estimated", "true", allocator);
            alpha.AddMember("min", this->alpha.min, allocator);
            alpha.AddMember("max", this->alpha.max, allocator);
            alpha.AddMember("phase", this->alpha.phase, allocator);

        } else {
            alpha.AddMember("estimated", "false", allocator);
            alpha.AddMember("min", this->alpha.min, allocator);
            alpha.AddMember("max", this->alpha.max, allocator);
            alpha.AddMember("phase", this->alpha.phase, allocator);
        }
        parameters.AddMember("alpha", alpha, allocator);

        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
            beta.AddMember("min", this->beta.min, allocator);
            beta.AddMember("max", this->beta.max, allocator);
            beta.AddMember("phase", this->beta.phase, allocator);

        } else {
            beta.AddMember("estimated", "false", allocator);
            beta.AddMember("min", this->beta.min, allocator);
            beta.AddMember("max", this->beta.max, allocator);
            beta.AddMember("phase", this->beta.phase, allocator);
        }
        parameters.AddMember("beta", beta, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kObjectType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        document.AddMember("recruitment", recruitment, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &rec, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "ricker", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value alpha(rapidjson::kObjectType);
        alpha.AddMember("value", this->alpha.value, allocator);
        if (this->alpha.estimated) {
            alpha.AddMember("estimated", "true", allocator);
            alpha.AddMember("min", this->alpha.min, allocator);
            alpha.AddMember("max", this->alpha.max, allocator);
            alpha.AddMember("phase", this->alpha.phase, allocator);

        } else {
            alpha.AddMember("estimated", "false", allocator);
            alpha.AddMember("min", this->alpha.min, allocator);
            alpha.AddMember("max", this->alpha.max, allocator);
            alpha.AddMember("phase", this->alpha.phase, allocator);
        }
        parameters.AddMember("alpha", alpha, allocator);

        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
            beta.AddMember("min", this->beta.min, allocator);
            beta.AddMember("max", this->beta.max, allocator);
            beta.AddMember("phase", this->beta.phase, allocator);

        } else {
            beta.AddMember("estimated", "false", allocator);
            beta.AddMember("min", this->beta.min, allocator);
            beta.AddMember("max", this->beta.max, allocator);
            beta.AddMember("phase", this->beta.phase, allocator);
        }
        parameters.AddMember("beta", beta, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kObjectType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        rec.PushBack(recruitment, allocator);

    }

    static std::map<int, RickerRecruitment*> initialized_models;
    typedef typename std::map<int, RickerRecruitment*>::iterator model_iterator;
};

std::map<int, RickerRecruitment*> RickerRecruitment::initialized_models;

class BevertonHoltRecruitment : public RecruitmentBase {
public:

    Rcpp::NumericVector deviations;
    double deviations_min = std::numeric_limits<double>::min();
    double deviations_max = std::numeric_limits<double>::max();
    bool estimate_deviations = true;
    int deviation_phase = 1;
    bool constrained_deviations = true;
    int id;
    Parameter R0;
    Parameter h;
    Parameter sigma_r;
    bool use_bias_correction = false;

    BevertonHoltRecruitment() {

        this->id = RecruitmentBase::id_g++;
        BevertonHoltRecruitment::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~BevertonHoltRecruitment() {
    }

    void SetDeviations(Rcpp::NumericVector values) {

        this->deviations = values;
    }

    double Evaluate(double SB0, double sb) {
        typedef typename mas::VariableTrait<double>::variable variable;
        variable ret;
        variable SB0_;
        variable sb_;
        mas::VariableTrait<double>::SetValue(SB0_, SB0);
        mas::VariableTrait<double>::SetValue(sb_, sb);
        atl::intrusive_ptr<mas::BevertonHolt<double> > rec =
                new mas::BevertonHolt<double>();
        mas::BevertonHolt<double> *r = rec.get();
        r->log_R0 = std::log(this->R0.value);
        r->R0 = this->R0.value;
        r->h = this->h.value;
        r->sigma_r = this->sigma_r.value;
        ret = r->Evaluate(SB0_, sb_);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        typedef typename mas::VariableTrait<double>::variable variable;

        atl::intrusive_ptr<mas::BevertonHolt<double> > rec =
                new mas::BevertonHolt<double>();
        mas::BevertonHolt<double> *r = rec.get();
        r->id = this->id;

        r->R0 = this->R0.value;
        r->log_R0 = std::log(this->R0.value);

        r->use_bias_correction = this->use_bias_correction;

        if (this->R0.estimated) {
            std::stringstream ss;
            ss << "log_R0_" << this->id;
            mas::VariableTrait<double>::SetName(r->log_R0, ss.str());
            if (this->R0.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->log_R0,
                        std::log(this->R0.min));
            }
            if (this->R0.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->log_R0,
                        std::log(this->R0.max));
            }
            r->Register(r->log_R0, this->R0.phase);
        }

        r->h = this->h.value;
        if (this->h.estimated) {
            std::stringstream ss;
            ss << "h" << this->id;
            mas::VariableTrait<double>::SetName(r->h, ss.str());

            if (this->h.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->h, h.min);
            }
            if (this->h.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->h, h.max);
            }
            r->Register(r->h, this->h.phase);
        }

        r->sigma_r = this->sigma_r.value;
        if (this->sigma_r.estimated) {
            std::stringstream ss;
            ss << "sigma_r" << this->id;
            mas::VariableTrait<double>::SetName(r->sigma_r, ss.str());
            if (this->sigma_r.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->sigma_r,
                        sigma_r.min);
            }
            if (this->sigma_r.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->sigma_r,
                        sigma_r.max);
            }
            r->Register(r->sigma_r, this->sigma_r.phase);
        }

        r->recruitment_deviations.resize(this->deviations.size());
        for (int i = 0; i < this->deviations.size(); i++) {
            r->recruitment_deviations[i] = variable(this->deviations[i]);
        }
        r->recruitment_deviations_constrained = this->constrained_deviations;

        if (this->estimate_deviations) {
            r->estimating_recruitment_deviations = true;
            for (int i = 0; i < info.nyears; i++) {
                std::stringstream ss;
                ss << "recruitment_deviations[" << i << "]_" << this->id;
                mas::VariableTrait<double>::SetName(
                        r->recruitment_deviations[i], ss.str());
                if (this->deviations_min
                        != std::numeric_limits<double>::min()) {
                    mas::VariableTrait<double>::SetMinBoundary(
                            r->recruitment_deviations[i], this->deviations_min);
                }
                if (this->deviations_max
                        != std::numeric_limits<double>::max()) {

                    mas::VariableTrait<double>::SetMaxBoundary(
                            r->recruitment_deviations[i], this->deviations_max);
                }
                r->Register(r->recruitment_deviations[i],
                        this->deviation_phase);

            }
        }

        info.recruitment_models[this->id] = rec;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::recruitment_model_iterator rit;
        rit = info.recruitment_models.find(this->id);
        if (rit != info.recruitment_models.end()) {

            atl::intrusive_ptr<mas::RecruitmentBase<double> > rec =
                    (*rit).second;
            mas::BevertonHolt<double> *r =
                    (mas::BevertonHolt<double>*) rec.get();
            this->R0.value = std::exp(r->log_R0.GetValue());
            this->h.value = r->h.GetValue();
            this->sigma_r.value = r->sigma_r.GetValue();
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "beverton_holt", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value h(rapidjson::kObjectType);
        h.AddMember("value", this->h.value, allocator);
        if (this->h.estimated) {
            h.AddMember("estimated", "true", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);

        } else {
            h.AddMember("estimated", "false", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);
        }
        parameters.AddMember("h", h, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        document.AddMember("recruitment", recruitment, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &rec, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "beverton_holt", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value h(rapidjson::kObjectType);
        h.AddMember("value", this->h.value, allocator);
        if (this->h.estimated) {
            h.AddMember("estimated", "true", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);

        } else {
            h.AddMember("estimated", "false", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);
        }
        parameters.AddMember("h", h, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        rec.PushBack(recruitment, allocator);
    }

    static std::map<int, BevertonHoltRecruitment*> initialized_models;
    typedef typename std::map<int, BevertonHoltRecruitment*>::iterator model_iterator;
};

std::map<int, BevertonHoltRecruitment*> BevertonHoltRecruitment::initialized_models;

class BevertonHoltRecruitmentAlt : public RecruitmentBase {
public:

    Rcpp::NumericVector deviations;
    double deviations_min = std::numeric_limits<double>::min();
    double deviations_max = std::numeric_limits<double>::max();
    bool estimate_deviations = true;
    int deviation_phase = 1;
    bool constrained_deviations = true;
    int id;
    Parameter R0;
    Parameter h;
    Parameter sigma_r;
    bool use_bias_correction = false;

    BevertonHoltRecruitmentAlt() {

        this->id = RecruitmentBase::id_g++;
        BevertonHoltRecruitmentAlt::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~BevertonHoltRecruitmentAlt() {
    }

    void SetDeviations(Rcpp::NumericVector values) {

        this->deviations = values;
    }

    double Evaluate(double SB0, double sb) {
        typedef typename mas::VariableTrait<double>::variable variable;
        variable ret;
        variable SB0_;
        variable sb_;
        mas::VariableTrait<double>::SetValue(SB0_, SB0);
        mas::VariableTrait<double>::SetValue(sb_, sb);
        atl::intrusive_ptr<mas::BevertonHolt<double> > rec =
                new mas::BevertonHolt<double>();
        mas::BevertonHolt<double> *r = rec.get();
        r->R0 = this->R0.value;
        r->R0 = this->R0.value;
        r->log_R0 = std::log(this->R0.value);
        r->h = this->h.value;
        r->sigma_r = this->sigma_r.value;
        ret = r->Evaluate(SB0_, sb_);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        typedef typename mas::VariableTrait<double>::variable variable;

        atl::intrusive_ptr<mas::BevertonHoltAlt<double> > rec =
                new mas::BevertonHoltAlt<double>();
        mas::BevertonHoltAlt<double> *r = rec.get();
        r->id = this->id;

        r->R0 = this->R0.value;
        r->log_R0 = std::log(this->R0.value);
        r->use_bias_correction = this->use_bias_correction;

        if (this->R0.estimated) {
            std::stringstream ss;
            ss << "log_R0_" << this->id;
            mas::VariableTrait<double>::SetName(r->log_R0, ss.str());
            if (this->R0.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->log_R0,
                        std::log(this->R0.min));
            }
            if (this->R0.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->log_R0,
                        std::log(this->R0.max));
            }
            r->Register(r->log_R0, this->R0.phase);
        }

        r->h = this->h.value;
        if (this->h.estimated) {
            std::stringstream ss;
            ss << "h" << this->id;
            mas::VariableTrait<double>::SetName(r->h, ss.str());

            if (this->h.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->h, h.min);
            }
            if (this->h.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->h, h.max);
            }
            r->Register(r->h, this->h.phase);
        }

        r->sigma_r = this->sigma_r.value;
        if (this->sigma_r.estimated) {
            std::stringstream ss;
            ss << "sigma_r" << this->id;
            mas::VariableTrait<double>::SetName(r->sigma_r, ss.str());
            if (this->sigma_r.min != std::numeric_limits<double>::min()) {
                mas::VariableTrait<double>::SetMinBoundary(r->sigma_r,
                        sigma_r.min);
            }
            if (this->sigma_r.max != std::numeric_limits<double>::max()) {
                mas::VariableTrait<double>::SetMaxBoundary(r->sigma_r,
                        sigma_r.max);
            }
            r->Register(r->sigma_r, this->sigma_r.phase);
        }

        r->recruitment_deviations.resize(this->deviations.size());
        for (int i = 0; i < this->deviations.size(); i++) {
            r->recruitment_deviations[i] = variable(this->deviations[i]);
        }
        r->recruitment_deviations_constrained = this->constrained_deviations;

        if (this->estimate_deviations) {
            r->estimating_recruitment_deviations = true;
            for (int i = 0; i < info.nyears; i++) {
                std::stringstream ss;
                ss << "recruitment_deviations[" << i << "]_" << this->id;
                mas::VariableTrait<double>::SetName(
                        r->recruitment_deviations[i], ss.str());
                if (this->deviations_min
                        != std::numeric_limits<double>::min()) {
                    mas::VariableTrait<double>::SetMinBoundary(
                            r->recruitment_deviations[i], this->deviations_min);
                }
                if (this->deviations_max
                        != std::numeric_limits<double>::max()) {

                    mas::VariableTrait<double>::SetMaxBoundary(
                            r->recruitment_deviations[i], this->deviations_max);
                }
                r->Register(r->recruitment_deviations[i],
                        this->deviation_phase);

            }
        }

        info.recruitment_models[this->id] = rec;
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::recruitment_model_iterator rit;
        rit = info.recruitment_models.find(this->id);
        if (rit != info.recruitment_models.end()) {

            atl::intrusive_ptr<mas::RecruitmentBase<double> > rec =
                    (*rit).second;
            mas::BevertonHoltAlt<double> *r =
                    (mas::BevertonHoltAlt<double>*) rec.get();
            this->R0.value = std::exp(r->log_R0.GetValue());
            this->h.value = r->h.GetValue();
            this->sigma_r.value = r->sigma_r.GetValue();
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "beverton_holt", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value h(rapidjson::kObjectType);
        h.AddMember("value", this->h.value, allocator);
        if (this->h.estimated) {
            h.AddMember("estimated", "true", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);

        } else {
            h.AddMember("estimated", "false", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);
        }
        parameters.AddMember("h", h, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        document.AddMember("recruitment", recruitment, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &rec, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value recruitment(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        recruitment.AddMember("id", this->id, allocator);
        recruitment.AddMember("model", "beverton_holt", allocator);

        rapidjson::Value R0(rapidjson::kObjectType);
        R0.AddMember("value", this->R0.value, allocator);
        if (this->R0.estimated) {
            R0.AddMember("estimated", "true", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);

        } else {
            R0.AddMember("estimated", "false", allocator);
            R0.AddMember("min", this->R0.min, allocator);
            R0.AddMember("max", this->R0.max, allocator);
            R0.AddMember("phase", this->R0.phase, allocator);
        }
        parameters.AddMember("R0", R0, allocator);

        rapidjson::Value h(rapidjson::kObjectType);
        h.AddMember("value", this->h.value, allocator);
        if (this->h.estimated) {
            h.AddMember("estimated", "true", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);

        } else {
            h.AddMember("estimated", "false", allocator);
            h.AddMember("min", this->h.min, allocator);
            h.AddMember("max", this->h.max, allocator);
            h.AddMember("phase", this->h.phase, allocator);
        }
        parameters.AddMember("h", h, allocator);

        rapidjson::Value sigma_r(rapidjson::kObjectType);
        sigma_r.AddMember("value", this->sigma_r.value, allocator);
        if (this->sigma_r.estimated) {
            sigma_r.AddMember("estimated", "true", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);

        } else {
            sigma_r.AddMember("estimated", "false", allocator);
            sigma_r.AddMember("min", this->sigma_r.min, allocator);
            sigma_r.AddMember("max", this->sigma_r.max, allocator);
            sigma_r.AddMember("phase", this->sigma_r.phase, allocator);
        }
        parameters.AddMember("sigma_r", sigma_r, allocator);

        rapidjson::Value rec_devs(rapidjson::kObjectType);
        rapidjson::Value vals(rapidjson::kArrayType);
        for (int i = 0; i < this->deviations.size(); i++) {
            vals.PushBack(this->deviations[i], allocator);
        }
        if (this->constrained_deviations) {
            rec_devs.AddMember("constrained", "true", allocator);
        } else {
            rec_devs.AddMember("constrained", "false", allocator);
        }
        rec_devs.AddMember("random_effect", "false", allocator);
        if (this->estimate_deviations) {
            rec_devs.AddMember("estimated", "true", allocator);
        } else {
            rec_devs.AddMember("estimated", "false", allocator);
        }
        rec_devs.AddMember("min", this->deviations_min, allocator);
        rec_devs.AddMember("max", this->deviations_max, allocator);
        rec_devs.AddMember("phase", this->deviation_phase, allocator);
        rec_devs.AddMember("values", vals, allocator);
        parameters.AddMember("recruitment_deviations", rec_devs, allocator);
        recruitment.AddMember("parameters", parameters, allocator);
        rec.PushBack(recruitment, allocator);
    }

    static std::map<int, BevertonHoltRecruitmentAlt*> initialized_models;
    typedef typename std::map<int, BevertonHoltRecruitmentAlt*>::iterator model_iterator;
};

std::map<int, BevertonHoltRecruitmentAlt*> BevertonHoltRecruitmentAlt::initialized_models;

/**
 * Growth
 */
class GrowthBase : public MASSubModel {
protected:

public:
    static int id_g;

    virtual ~GrowthBase() {
    }
};
int GrowthBase::id_g = 1;

class VonBertalanffy : public GrowthBase {
public:

    Rcpp::NumericVector weight_at_season_start_females;
    Rcpp::NumericVector weight_at_spawning_females;
    Rcpp::NumericVector weight_at_catch_females;
    Rcpp::NumericVector weight_at_survey_females;
    Rcpp::NumericVector weight_at_season_start_males;
    Rcpp::NumericVector weight_at_spawning_males;
    Rcpp::NumericVector weight_at_catch_males;
    Rcpp::NumericVector weight_at_survey_males;
    bool has_emprical_weight = false;
    Parameter a_min;
    Parameter a_max;
    Parameter alpha_f = 0.000025;
    Parameter alpha_m = 0.000025;
    Parameter beta_f = 3.0;
    Parameter beta_m = 3.0;
    Parameter k;
    Parameter l_inf;
    int id;

    VonBertalanffy() {

        this->id = GrowthBase::id_g++;
        VonBertalanffy::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~VonBertalanffy() {
    }

    void SetUndifferentiatedWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->weight_at_season_start_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->weight_at_spawning_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->weight_at_catch_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->weight_at_survey_females = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->has_emprical_weight = true;
    }

    double Evaluate(const double &age, const int &sex) {
        atl::intrusive_ptr<mas::VonBertalanffy<double> > vb =
                new mas::VonBertalanffy<double>();
        mas::VonBertalanffy<double> *g = vb.get();
        mas::VariableTrait<double>::variable a;
        mas::VariableTrait<double>::SetValue(a, age);
        mas::VariableTrait<double>::variable ret = g->Evaluate(a, sex);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::VonBertalanffy<double> > vb =
                new mas::VonBertalanffy<double>();
        mas::VonBertalanffy<double> *g = vb.get();
        g->id = this->id;
        std::stringstream ss;

        ss << "amin_" << this->id;
        this->InitializeParameter(g, g->a_min, this->a_min, ss.str());
        ss.str("");
        ss << "amax_" << this->id;
        this->InitializeParameter(g, g->a_max, this->a_max, ss.str());
        ss.str("");
        ss << "alpha_f_" << this->id;
        this->InitializeParameter(g, g->alpha_f, this->alpha_f, ss.str());
        ss.str("");
        ss << "alpha_m_" << this->id;
        this->InitializeParameter(g, g->alpha_m, this->alpha_m, ss.str());
        ss.str("");
        ss << "beta_f_" << this->id;
        this->InitializeParameter(g, g->beta_f, this->beta_f, ss.str());
        ss.str("");
        ss << "beta_m_" << this->id;
        this->InitializeParameter(g, g->beta_m, this->beta_m, ss.str());
        ss.str("");
        ss << "k_" << this->id;
        this->InitializeParameter(g, g->k, this->k, ss.str());
        ss.str("");
        ss << "l_inf_" << this->id;
        this->InitializeParameter(g, g->l_inf, this->l_inf, ss.str());

        size_t data_length = info.ages.size() * info.nyears * info.nseasons;

        atl::intrusive_ptr<mas::WeightFunctorBase<double> > weight_functor;
        if (this->has_emprical_weight == true) {

            weight_functor = new mas::EmpiricalWeightFunctor<double>();
            g->weight_functor = weight_functor;

            mas::EmpiricalWeightFunctor<double> *eg =
                    (mas::EmpiricalWeightFunctor<double>*) weight_functor.get();

            /**
             * Catch weight
             */
            if (this->weight_at_catch_females.size()) {
                if (this->weight_at_catch_females.size() != data_length) {
                    std::cout
                            << "weight_at_catch_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_catch_females[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::CATCH_MEAN_WEIGHT_AT_AGE][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_catch_males.size()) {
                if (this->weight_at_catch_males.size() != data_length) {
                    std::cout
                            << "weight_at_catch_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_catch_males[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::CATCH_MEAN_WEIGHT_AT_AGE][mas::MALE] =
                        eds;
            }
            /*
             * Survey weight
             */
            if (this->weight_at_survey_females.size()) {
                if (this->weight_at_survey_females.size() != data_length) {
                    std::cout
                            << "weight_at_survey_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_survey_females[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::SURVEY_MEAN_WEIGHT_AT_AGE][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_survey_males.size()) {
                if (this->weight_at_survey_males.size() != data_length) {
                    std::cout
                            << "weight_at_survey_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_survey_males[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::SURVEY_MEAN_WEIGHT_AT_AGE][mas::MALE] =
                        eds;
            }
            /*
             * Spawning weight
             */
            if (this->weight_at_spawning_females.size()) {
                if (this->weight_at_spawning_females.size() != data_length) {
                    std::cout
                            << "weight_at_spawning_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_spawning_females[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SPAWNING][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_spawning_males.size()) {
                if (this->weight_at_spawning_males.size() != data_length) {
                    std::cout
                            << "weight_at_spawning_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_spawning_males[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SPAWNING][mas::MALE] =
                        eds;
            }
            /*
             * Season start weight
             */
            if (this->weight_at_season_start_females.size()) {
                if (this->weight_at_season_start_females.size()
                        != data_length) {
                    std::cout
                            << "weight_at_season_start_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_season_start_females[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SEASON_START][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_season_start_males.size()) {
                if (this->weight_at_season_start_males.size() != data_length) {
                    std::cout
                            << "weight_at_season_start_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_season_start_males[i];
                }
                eds.empirical_data_at_age = data;
                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SEASON_START][mas::MALE] =
                        eds;
            }
        } else {

            //instantiate default here

            weight_functor = new mas::DefaultWeightFunctor<double>(g->alpha_f,
                    g->alpha_m, g->beta_f, g->beta_f);
            g->weight_functor = weight_functor;
        }

        info.growth_models[this->id] = vb;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value growth(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        growth.AddMember("id", this->id, allocator);
        growth.AddMember("model", "von_bertalanffy", allocator);

        //amin
        rapidjson::Value amin(rapidjson::kObjectType);
        amin.AddMember("value", this->a_min.value, allocator);
        if (this->a_min.estimated) {
            amin.AddMember("estimated", "true", allocator);
        } else {
            amin.AddMember("estimated", "false", allocator);
        }
        amin.AddMember("min", this->a_min.min, allocator);
        amin.AddMember("min", this->a_min.max, allocator);
        amin.AddMember("phase", this->a_min.phase, allocator);
        parameters.AddMember("amin", amin, allocator);

        //amax
        rapidjson::Value amax(rapidjson::kObjectType);
        amax.AddMember("value", this->a_max.value, allocator);
        if (this->a_max.estimated) {
            amax.AddMember("estimated", "true", allocator);
        } else {
            amax.AddMember("estimated", "false", allocator);
        }
        amax.AddMember("min", this->a_max.min, allocator);
        amax.AddMember("min", this->a_max.max, allocator);
        amax.AddMember("phase", this->a_max.phase, allocator);
        parameters.AddMember("amax", amax, allocator);

        //k
        rapidjson::Value k(rapidjson::kObjectType);
        k.AddMember("value", this->k.value, allocator);
        if (this->k.estimated) {
            k.AddMember("estimated", "true", allocator);
        } else {
            k.AddMember("estimated", "false", allocator);
        }
        k.AddMember("min", this->k.min, allocator);
        k.AddMember("min", this->k.max, allocator);
        k.AddMember("phase", this->k.phase, allocator);
        parameters.AddMember("k", k, allocator);

        //linf
        rapidjson::Value linf(rapidjson::kObjectType);
        linf.AddMember("value", this->l_inf.value, allocator);
        if (this->l_inf.estimated) {
            linf.AddMember("estimated", "true", allocator);
        } else {
            linf.AddMember("estimated", "false", allocator);
        }
        linf.AddMember("min", this->l_inf.min, allocator);
        linf.AddMember("min", this->l_inf.max, allocator);
        linf.AddMember("phase", this->l_inf.phase, allocator);
        parameters.AddMember("linf", linf, allocator);

        //alpha f
        rapidjson::Value alpha_f(rapidjson::kObjectType);
        alpha_f.AddMember("value", this->alpha_f.value, allocator);
        if (this->alpha_f.estimated) {
            alpha_f.AddMember("estimated", "true", allocator);
        } else {
            alpha_f.AddMember("estimated", "false", allocator);
        }
        alpha_f.AddMember("min", this->alpha_f.min, allocator);
        alpha_f.AddMember("min", this->alpha_f.max, allocator);
        alpha_f.AddMember("phase", this->alpha_f.phase, allocator);
        parameters.AddMember("alpha_f", alpha_f, allocator);

        //alpha m
        rapidjson::Value alpha_m(rapidjson::kObjectType);
        alpha_m.AddMember("value", this->alpha_m.value, allocator);
        if (this->alpha_m.estimated) {
            alpha_m.AddMember("estimated", "true", allocator);
        } else {
            alpha_m.AddMember("estimated", "false", allocator);
        }
        alpha_m.AddMember("min", this->alpha_m.min, allocator);
        alpha_m.AddMember("min", this->alpha_m.max, allocator);
        alpha_m.AddMember("phase", this->alpha_m.phase, allocator);
        parameters.AddMember("alpha_m", alpha_m, allocator);

        //beta f
        rapidjson::Value beta_f(rapidjson::kObjectType);
        beta_f.AddMember("value", this->beta_f.value, allocator);
        if (this->beta_f.estimated) {
            beta_f.AddMember("estimated", "true", allocator);
        } else {
            beta_f.AddMember("estimated", "false", allocator);
        }
        beta_f.AddMember("min", this->beta_f.min, allocator);
        beta_f.AddMember("min", this->beta_f.max, allocator);
        beta_f.AddMember("phase", this->beta_f.phase, allocator);
        parameters.AddMember("beta_f", beta_f, allocator);

        //beta m
        rapidjson::Value beta_m(rapidjson::kObjectType);
        beta_m.AddMember("value", this->beta_m.value, allocator);
        if (this->a_max.estimated) {
            beta_m.AddMember("estimated", "true", allocator);
        } else {
            beta_m.AddMember("estimated", "false", allocator);
        }
        beta_m.AddMember("min", this->beta_m.min, allocator);
        beta_m.AddMember("min", this->beta_m.max, allocator);
        beta_m.AddMember("phase", this->beta_m.phase, allocator);
        parameters.AddMember("beta_m", beta_m, allocator);

        growth.AddMember("parameters", parameters, allocator);

        //empirical weight at age
        if (this->has_emprical_weight) {
            rapidjson::Value empirical_weight_at_age(rapidjson::kArrayType);

            if (this->weight_at_catch_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_catch_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }
            growth.AddMember("empirical_weight_at_age", empirical_weight_at_age,
                    allocator);

        }

        document.AddMember("growth", growth, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &grth, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value growth(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        growth.AddMember("id", this->id, allocator);
        growth.AddMember("model", "von_bertalanffy", allocator);

        //amin
        rapidjson::Value amin(rapidjson::kObjectType);
        amin.AddMember("value", this->a_min.value, allocator);
        if (this->a_min.estimated) {
            amin.AddMember("estimated", "true", allocator);
        } else {
            amin.AddMember("estimated", "false", allocator);
        }
        amin.AddMember("min", this->a_min.min, allocator);
        amin.AddMember("min", this->a_min.max, allocator);
        amin.AddMember("phase", this->a_min.phase, allocator);
        parameters.AddMember("amin", amin, allocator);

        //amax
        rapidjson::Value amax(rapidjson::kObjectType);
        amax.AddMember("value", this->a_max.value, allocator);
        if (this->a_max.estimated) {
            amax.AddMember("estimated", "true", allocator);
        } else {
            amax.AddMember("estimated", "false", allocator);
        }
        amax.AddMember("min", this->a_max.min, allocator);
        amax.AddMember("min", this->a_max.max, allocator);
        amax.AddMember("phase", this->a_max.phase, allocator);
        parameters.AddMember("amax", amax, allocator);

        //k
        rapidjson::Value k(rapidjson::kObjectType);
        k.AddMember("value", this->k.value, allocator);
        if (this->k.estimated) {
            k.AddMember("estimated", "true", allocator);
        } else {
            k.AddMember("estimated", "false", allocator);
        }
        k.AddMember("min", this->k.min, allocator);
        k.AddMember("min", this->k.max, allocator);
        k.AddMember("phase", this->k.phase, allocator);
        parameters.AddMember("k", k, allocator);

        //linf
        rapidjson::Value linf(rapidjson::kObjectType);
        linf.AddMember("value", this->l_inf.value, allocator);
        if (this->l_inf.estimated) {
            linf.AddMember("estimated", "true", allocator);
        } else {
            linf.AddMember("estimated", "false", allocator);
        }
        linf.AddMember("min", this->l_inf.min, allocator);
        linf.AddMember("min", this->l_inf.max, allocator);
        linf.AddMember("phase", this->l_inf.phase, allocator);
        parameters.AddMember("linf", linf, allocator);

        //alpha f
        rapidjson::Value alpha_f(rapidjson::kObjectType);
        alpha_f.AddMember("value", this->alpha_f.value, allocator);
        if (this->alpha_f.estimated) {
            alpha_f.AddMember("estimated", "true", allocator);
        } else {
            alpha_f.AddMember("estimated", "false", allocator);
        }
        alpha_f.AddMember("min", this->alpha_f.min, allocator);
        alpha_f.AddMember("min", this->alpha_f.max, allocator);
        alpha_f.AddMember("phase", this->alpha_f.phase, allocator);
        parameters.AddMember("alpha_f", alpha_f, allocator);

        //alpha m
        rapidjson::Value alpha_m(rapidjson::kObjectType);
        alpha_m.AddMember("value", this->alpha_m.value, allocator);
        if (this->alpha_m.estimated) {
            alpha_m.AddMember("estimated", "true", allocator);
        } else {
            alpha_m.AddMember("estimated", "false", allocator);
        }
        alpha_m.AddMember("min", this->alpha_m.min, allocator);
        alpha_m.AddMember("min", this->alpha_m.max, allocator);
        alpha_m.AddMember("phase", this->alpha_m.phase, allocator);
        parameters.AddMember("alpha_m", alpha_m, allocator);

        //beta f
        rapidjson::Value beta_f(rapidjson::kObjectType);
        beta_f.AddMember("value", this->beta_f.value, allocator);
        if (this->beta_f.estimated) {
            beta_f.AddMember("estimated", "true", allocator);
        } else {
            beta_f.AddMember("estimated", "false", allocator);
        }
        beta_f.AddMember("min", this->beta_f.min, allocator);
        beta_f.AddMember("min", this->beta_f.max, allocator);
        beta_f.AddMember("phase", this->beta_f.phase, allocator);
        parameters.AddMember("beta_f", beta_f, allocator);

        //beta m
        rapidjson::Value beta_m(rapidjson::kObjectType);
        beta_m.AddMember("value", this->beta_m.value, allocator);
        if (this->a_max.estimated) {
            beta_m.AddMember("estimated", "true", allocator);
        } else {
            beta_m.AddMember("estimated", "false", allocator);
        }
        beta_m.AddMember("min", this->beta_m.min, allocator);
        beta_m.AddMember("min", this->beta_m.max, allocator);
        beta_m.AddMember("phase", this->beta_m.phase, allocator);
        parameters.AddMember("beta_m", beta_m, allocator);

        growth.AddMember("parameters", parameters, allocator);

        //empirical weight at age
        if (this->has_emprical_weight) {
            rapidjson::Value empirical_weight_at_age(rapidjson::kArrayType);

            if (this->weight_at_catch_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_catch_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }
            growth.AddMember("empirical_weight_at_age", empirical_weight_at_age,
                    allocator);

        }

        grth.PushBack(growth, allocator);
    }

    static std::map<int, VonBertalanffy*> initialized_models;
    typedef typename std::map<int, VonBertalanffy*>::iterator model_iterator;
};

std::map<int, VonBertalanffy*> VonBertalanffy::initialized_models;

class VonBertalanffyModified : public GrowthBase {
public:

    Rcpp::NumericVector weight_at_season_start_females;
    Rcpp::NumericVector weight_at_spawning_females;
    Rcpp::NumericVector weight_at_catch_females;
    Rcpp::NumericVector weight_at_survey_females;
    Rcpp::NumericVector weight_at_season_start_males;
    Rcpp::NumericVector weight_at_spawning_males;
    Rcpp::NumericVector weight_at_catch_males;
    Rcpp::NumericVector weight_at_survey_males;
    bool has_emprical_weight = false;
    Parameter a_min;
    Parameter a_max;
    Parameter alpha_f = 0.000025;
    Parameter alpha_m = 0.000025;
    Parameter beta_f = 3.0;
    Parameter beta_m = 3.0;
    Parameter lmin;
    Parameter lmax;
    Parameter l_inf;
    Parameter c;
    int id;

    VonBertalanffyModified() {

        this->id = GrowthBase::id_g++;
        VonBertalanffyModified::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~VonBertalanffyModified() {
    }

    void SetUndifferentiatedWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->weight_at_season_start_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->weight_at_spawning_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->weight_at_catch_females = data;
        this->has_emprical_weight = true;
    }

    void SetUndifferentiatedSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->weight_at_survey_females = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->has_emprical_weight = true;
    }

    void SetFemaleSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleWeightAtSeasonStart(Rcpp::NumericVector data) {

        this->weight_at_season_start_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleWeightAtSpawning(Rcpp::NumericVector data) {

        this->weight_at_spawning_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleCatchWeight(Rcpp::NumericVector data) {

        this->weight_at_catch_males = data;
        this->has_emprical_weight = true;
    }

    void SetMaleSurveyWeight(Rcpp::NumericVector data) {

        this->weight_at_survey_males = data;
        this->has_emprical_weight = true;
    }

    double Evaluate(const double &age, const int &sex) {
        atl::intrusive_ptr<mas::VonBertalanffyModified<double> > vb =
                new mas::VonBertalanffyModified<double>();
        mas::VonBertalanffyModified<double> *g = vb.get();
        mas::VariableTrait<double>::variable a;
        mas::VariableTrait<double>::SetValue(a, age);
        mas::VariableTrait<double>::variable ret = g->Evaluate(a, sex);
        return mas::VariableTrait<double>::Value(ret);
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::VonBertalanffyModified<double> > vb =
                new mas::VonBertalanffyModified<double>();
        mas::VonBertalanffyModified<double> *g = vb.get();
        g->id = this->id;
        std::stringstream ss;

        ss << "amin_" << this->id;
        this->InitializeParameter(g, g->a_min, this->a_min, ss.str());
        ss.str("");
        ss << "amax_" << this->id;
        this->InitializeParameter(g, g->a_max, this->a_max, ss.str());
        ss.str("");
        ss << "alpha_f_" << this->id;
        this->InitializeParameter(g, g->alpha_f, this->alpha_f, ss.str());
        ss.str("");
        ss << "alpha_m_" << this->id;
        this->InitializeParameter(g, g->alpha_m, this->alpha_m, ss.str());
        ss.str("");
        ss << "beta_f_" << this->id;
        this->InitializeParameter(g, g->beta_f, this->beta_f, ss.str());
        ss.str("");
        ss << "beta_m_" << this->id;
        this->InitializeParameter(g, g->beta_m, this->beta_m, ss.str());
        ss.str("");
        ss << "lmin_" << this->id;
        this->InitializeParameter(g, g->lmin, this->lmin, ss.str());
        ss.str("");
        ss << "lmax_" << this->id;
        this->InitializeParameter(g, g->lmax, this->lmax, ss.str());
        ss.str("");
        ss << "c_" << this->id;
        this->InitializeParameter(g, g->c, this->c, ss.str());
        ss.str("");
        ss << "l_inf_" << this->id;
        this->InitializeParameter(g, g->l_inf, this->l_inf, ss.str());

        atl::intrusive_ptr<mas::WeightFunctorBase<double> > weight_functor;
        size_t data_length = info.ages.size() * info.nyears * info.nseasons;
        if (this->has_emprical_weight == true) {

            weight_functor = new mas::EmpiricalWeightFunctor<double>();
            //            g->weight_functor = weight_functor;

            mas::EmpiricalWeightFunctor<double> *eg =
                    (mas::EmpiricalWeightFunctor<double>*) weight_functor.get();

            /**
             * Catch weight
             */
            if (this->weight_at_catch_females.size()) {
                if (this->weight_at_catch_females.size() != data_length) {
                    std::cout
                            << "weight_at_catch_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->type = mas::CATCH_MEAN_WEIGHT_AT_AGE;
                //                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data.push_back(this->weight_at_catch_females[i]);
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::CATCH_MEAN_WEIGHT_AT_AGE][mas::FEMALE] =
                        eds;

            }

            if (this->weight_at_catch_males.size()) {
                if (this->weight_at_catch_males.size() != data_length) {
                    std::cout
                            << "weight_at_catch_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_catch_males[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::CATCH_MEAN_WEIGHT_AT_AGE][mas::MALE] =
                        eds;
            }
            /*
             * Survey weight
             */
            if (this->weight_at_survey_females.size()) {
                if (this->weight_at_survey_females.size() != data_length) {
                    std::cout
                            << "weight_at_survey_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_survey_females[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::SURVEY_MEAN_WEIGHT_AT_AGE][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_survey_males.size()) {
                if (this->weight_at_survey_males.size() != data_length) {
                    std::cout
                            << "weight_at_survey_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_survey_males[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::SURVEY_MEAN_WEIGHT_AT_AGE][mas::MALE] =
                        eds;
            }
            /*
             * Spawning weight
             */
            if (this->weight_at_spawning_females.size()) {
                if (this->weight_at_spawning_females.size() != data_length) {
                    std::cout
                            << "weight_at_spawning_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_spawning_females[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SPAWNING][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_spawning_males.size()) {
                if (this->weight_at_spawning_males.size() != data_length) {
                    std::cout
                            << "weight_at_spawning_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_spawning_males[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SPAWNING][mas::MALE] =
                        eds;
            }
            /*
             * Season start weight
             */
            if (this->weight_at_season_start_females.size()) {
                if (this->weight_at_season_start_females.size()
                        != data_length) {
                    std::cout
                            << "weight_at_season_start_females vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_season_start_females[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SEASON_START][mas::FEMALE] =
                        eds;
            }

            if (this->weight_at_season_start_males.size()) {
                if (this->weight_at_season_start_males.size() != data_length) {
                    std::cout
                            << "weight_at_season_start_males vector not equal to (nyears*nseasons*nages)\n";
                    info.valid_configuration = false;
                    return;
                }
                size_t ages = info.ages.size();
                size_t years = info.nyears;
                size_t seasons = info.nseasons;

                mas::EmpricalDataStructure<double> eds;
                atl::intrusive_ptr<mas::DataObject<double> > data =
                        new mas::DataObject<double>();
                mas::DataObject<double> *d = data.get();
                d->imax = years;
                d->jmax = seasons;
                d->kmax = ages;
                d->dimensions = 3;
                d->data.resize(data_length);
                for (int i = 0; i < data_length; i++) {
                    d->data[i] = this->weight_at_season_start_males[i];
                }
                eds.empirical_data_at_age = data;

                mas::GrowthBase<double>::Do3DInterpolation(
                        eds.empirical_data_at_age,
                        eds.interpolated_data_at_age);
                eg->weight_at_age_data[mas::MEAN_WEIGHT_AT_AGE_SEASON_START][mas::MALE] =
                        eds;
            }
        } else {

            //instantiate default here

            weight_functor = new mas::DefaultWeightFunctor<double>(g->alpha_f,
                    g->alpha_m, g->beta_f, g->beta_f);
        }
        vb->weight_functor = weight_functor;
        info.growth_models[this->id] = vb;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value growth(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        growth.AddMember("id", this->id, allocator);
        growth.AddMember("model", "von_bertalanffy_modified", allocator);
        /*
         Parameter a_min;
         Parameter a_max;
         Parameter alpha_f = 0.000025;
         Parameter alpha_m = 0.000025;
         Parameter beta_f = 3.0;
         Parameter beta_m = 3.0;
         Parameter lmin;
         Parameter lmax;
         Parameter l_inf;
         Parameter c;
         */

        //amin
        rapidjson::Value amin(rapidjson::kObjectType);
        amin.AddMember("value", this->a_min.value, allocator);
        if (this->a_min.estimated) {
            amin.AddMember("estimated", "true", allocator);
        } else {
            amin.AddMember("estimated", "false", allocator);
        }
        amin.AddMember("min", this->a_min.min, allocator);
        amin.AddMember("min", this->a_min.max, allocator);
        amin.AddMember("phase", this->a_min.phase, allocator);
        parameters.AddMember("amin", amin, allocator);

        //amax
        rapidjson::Value amax(rapidjson::kObjectType);
        amax.AddMember("value", this->a_max.value, allocator);
        if (this->a_max.estimated) {
            amax.AddMember("estimated", "true", allocator);
        } else {
            amax.AddMember("estimated", "false", allocator);
        }
        amax.AddMember("min", this->a_max.min, allocator);
        amax.AddMember("min", this->a_max.max, allocator);
        amax.AddMember("phase", this->a_max.phase, allocator);
        parameters.AddMember("amax", amax, allocator);

        //lmin
        rapidjson::Value lmin(rapidjson::kObjectType);
        lmin.AddMember("value", this->lmin.value, allocator);
        if (this->lmin.estimated) {
            lmin.AddMember("estimated", "true", allocator);
        } else {
            lmin.AddMember("estimated", "false", allocator);
        }
        lmin.AddMember("min", this->lmin.min, allocator);
        lmin.AddMember("min", this->lmin.max, allocator);
        lmin.AddMember("phase", this->lmin.phase, allocator);
        parameters.AddMember("lmin", lmin, allocator);

        //lmax
        rapidjson::Value lmax(rapidjson::kObjectType);
        lmax.AddMember("value", this->lmax.value, allocator);
        if (this->lmax.estimated) {
            lmax.AddMember("estimated", "true", allocator);
        } else {
            lmax.AddMember("estimated", "false", allocator);
        }
        lmax.AddMember("min", this->lmax.min, allocator);
        lmax.AddMember("min", this->lmax.max, allocator);
        lmax.AddMember("phase", this->lmax.phase, allocator);
        parameters.AddMember("lmax", lmax, allocator);

        //c
        rapidjson::Value c(rapidjson::kObjectType);
        c.AddMember("value", this->c.value, allocator);
        if (this->c.estimated) {
            c.AddMember("estimated", "true", allocator);
        } else {
            c.AddMember("estimated", "false", allocator);
        }
        c.AddMember("min", this->c.min, allocator);
        c.AddMember("min", this->c.max, allocator);
        c.AddMember("phase", this->c.phase, allocator);
        parameters.AddMember("c", c, allocator);

        //alpha f
        rapidjson::Value alpha_f(rapidjson::kObjectType);
        alpha_f.AddMember("value", this->alpha_f.value, allocator);
        if (this->alpha_f.estimated) {
            alpha_f.AddMember("estimated", "true", allocator);
        } else {
            alpha_f.AddMember("estimated", "false", allocator);
        }
        alpha_f.AddMember("min", this->alpha_f.min, allocator);
        alpha_f.AddMember("min", this->alpha_f.max, allocator);
        alpha_f.AddMember("phase", this->alpha_f.phase, allocator);
        parameters.AddMember("alpha_f", alpha_f, allocator);

        //alpha m
        rapidjson::Value alpha_m(rapidjson::kObjectType);
        alpha_m.AddMember("value", this->alpha_m.value, allocator);
        if (this->alpha_m.estimated) {
            alpha_m.AddMember("estimated", "true", allocator);
        } else {
            alpha_m.AddMember("estimated", "false", allocator);
        }
        alpha_m.AddMember("min", this->alpha_m.min, allocator);
        alpha_m.AddMember("min", this->alpha_m.max, allocator);
        alpha_m.AddMember("phase", this->alpha_m.phase, allocator);
        parameters.AddMember("alpha_m", alpha_m, allocator);

        //beta f
        rapidjson::Value beta_f(rapidjson::kObjectType);
        beta_f.AddMember("value", this->beta_f.value, allocator);
        if (this->beta_f.estimated) {
            beta_f.AddMember("estimated", "true", allocator);
        } else {
            beta_f.AddMember("estimated", "false", allocator);
        }
        beta_f.AddMember("min", this->beta_f.min, allocator);
        beta_f.AddMember("min", this->beta_f.max, allocator);
        beta_f.AddMember("phase", this->beta_f.phase, allocator);
        parameters.AddMember("beta_f", beta_f, allocator);

        //beta m
        rapidjson::Value beta_m(rapidjson::kObjectType);
        beta_m.AddMember("value", this->beta_m.value, allocator);
        if (this->a_max.estimated) {
            beta_m.AddMember("estimated", "true", allocator);
        } else {
            beta_m.AddMember("estimated", "false", allocator);
        }
        beta_m.AddMember("min", this->beta_m.min, allocator);
        beta_m.AddMember("min", this->beta_m.max, allocator);
        beta_m.AddMember("phase", this->beta_m.phase, allocator);
        parameters.AddMember("beta_m", beta_m, allocator);

        //linf
        rapidjson::Value linf(rapidjson::kObjectType);
        linf.AddMember("value", this->l_inf.value, allocator);
        if (this->l_inf.estimated) {
            linf.AddMember("estimated", "true", allocator);
        } else {
            linf.AddMember("estimated", "false", allocator);
        }
        linf.AddMember("min", this->l_inf.min, allocator);
        linf.AddMember("min", this->l_inf.max, allocator);
        linf.AddMember("phase", this->l_inf.phase, allocator);
        parameters.AddMember("linf", linf, allocator);

        growth.AddMember("parameters", parameters, allocator);

        //empirical weight at age
        if (this->has_emprical_weight) {
            rapidjson::Value empirical_weight_at_age(rapidjson::kArrayType);

            if (this->weight_at_catch_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_catch_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }
            growth.AddMember("empirical_weight_at_age", empirical_weight_at_age,
                    allocator);

        }
        document.AddMember("growth", growth, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &grth, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value growth(rapidjson::kObjectType);
        rapidjson::Value parameters(rapidjson::kObjectType);

        growth.AddMember("id", this->id, allocator);
        growth.AddMember("model", "von_bertalanffy_modified", allocator);
        /*
         Parameter a_min;
         Parameter a_max;
         Parameter alpha_f = 0.000025;
         Parameter alpha_m = 0.000025;
         Parameter beta_f = 3.0;
         Parameter beta_m = 3.0;
         Parameter lmin;
         Parameter lmax;
         Parameter l_inf;
         Parameter c;
         */

        //amin
        rapidjson::Value amin(rapidjson::kObjectType);
        amin.AddMember("value", this->a_min.value, allocator);
        if (this->a_min.estimated) {
            amin.AddMember("estimated", "true", allocator);
        } else {
            amin.AddMember("estimated", "false", allocator);
        }
        amin.AddMember("min", this->a_min.min, allocator);
        amin.AddMember("min", this->a_min.max, allocator);
        amin.AddMember("phase", this->a_min.phase, allocator);
        parameters.AddMember("amin", amin, allocator);

        //amax
        rapidjson::Value amax(rapidjson::kObjectType);
        amax.AddMember("value", this->a_max.value, allocator);
        if (this->a_max.estimated) {
            amax.AddMember("estimated", "true", allocator);
        } else {
            amax.AddMember("estimated", "false", allocator);
        }
        amax.AddMember("min", this->a_max.min, allocator);
        amax.AddMember("min", this->a_max.max, allocator);
        amax.AddMember("phase", this->a_max.phase, allocator);
        parameters.AddMember("amax", amax, allocator);

        //lmin
        rapidjson::Value lmin(rapidjson::kObjectType);
        lmin.AddMember("value", this->lmin.value, allocator);
        if (this->lmin.estimated) {
            lmin.AddMember("estimated", "true", allocator);
        } else {
            lmin.AddMember("estimated", "false", allocator);
        }
        lmin.AddMember("min", this->lmin.min, allocator);
        lmin.AddMember("min", this->lmin.max, allocator);
        lmin.AddMember("phase", this->lmin.phase, allocator);
        parameters.AddMember("lmin", lmin, allocator);

        //lmax
        rapidjson::Value lmax(rapidjson::kObjectType);
        lmax.AddMember("value", this->lmax.value, allocator);
        if (this->lmax.estimated) {
            lmax.AddMember("estimated", "true", allocator);
        } else {
            lmax.AddMember("estimated", "false", allocator);
        }
        lmax.AddMember("min", this->lmax.min, allocator);
        lmax.AddMember("min", this->lmax.max, allocator);
        lmax.AddMember("phase", this->lmax.phase, allocator);
        parameters.AddMember("lmax", lmax, allocator);

        //c
        rapidjson::Value c(rapidjson::kObjectType);
        c.AddMember("value", this->c.value, allocator);
        if (this->c.estimated) {
            c.AddMember("estimated", "true", allocator);
        } else {
            c.AddMember("estimated", "false", allocator);
        }
        c.AddMember("min", this->c.min, allocator);
        c.AddMember("min", this->c.max, allocator);
        c.AddMember("phase", this->c.phase, allocator);
        parameters.AddMember("c", c, allocator);

        //alpha f
        rapidjson::Value alpha_f(rapidjson::kObjectType);
        alpha_f.AddMember("value", this->alpha_f.value, allocator);
        if (this->alpha_f.estimated) {
            alpha_f.AddMember("estimated", "true", allocator);
        } else {
            alpha_f.AddMember("estimated", "false", allocator);
        }
        alpha_f.AddMember("min", this->alpha_f.min, allocator);
        alpha_f.AddMember("min", this->alpha_f.max, allocator);
        alpha_f.AddMember("phase", this->alpha_f.phase, allocator);
        parameters.AddMember("alpha_f", alpha_f, allocator);

        //alpha m
        rapidjson::Value alpha_m(rapidjson::kObjectType);
        alpha_m.AddMember("value", this->alpha_m.value, allocator);
        if (this->alpha_m.estimated) {
            alpha_m.AddMember("estimated", "true", allocator);
        } else {
            alpha_m.AddMember("estimated", "false", allocator);
        }
        alpha_m.AddMember("min", this->alpha_m.min, allocator);
        alpha_m.AddMember("min", this->alpha_m.max, allocator);
        alpha_m.AddMember("phase", this->alpha_m.phase, allocator);
        parameters.AddMember("alpha_m", alpha_m, allocator);

        //beta f
        rapidjson::Value beta_f(rapidjson::kObjectType);
        beta_f.AddMember("value", this->beta_f.value, allocator);
        if (this->beta_f.estimated) {
            beta_f.AddMember("estimated", "true", allocator);
        } else {
            beta_f.AddMember("estimated", "false", allocator);
        }
        beta_f.AddMember("min", this->beta_f.min, allocator);
        beta_f.AddMember("min", this->beta_f.max, allocator);
        beta_f.AddMember("phase", this->beta_f.phase, allocator);
        parameters.AddMember("beta_f", beta_f, allocator);

        //beta m
        rapidjson::Value beta_m(rapidjson::kObjectType);
        beta_m.AddMember("value", this->beta_m.value, allocator);
        if (this->a_max.estimated) {
            beta_m.AddMember("estimated", "true", allocator);
        } else {
            beta_m.AddMember("estimated", "false", allocator);
        }
        beta_m.AddMember("min", this->beta_m.min, allocator);
        beta_m.AddMember("min", this->beta_m.max, allocator);
        beta_m.AddMember("phase", this->beta_m.phase, allocator);
        parameters.AddMember("beta_m", beta_m, allocator);

        //linf
        rapidjson::Value linf(rapidjson::kObjectType);
        linf.AddMember("value", this->l_inf.value, allocator);
        if (this->l_inf.estimated) {
            linf.AddMember("estimated", "true", allocator);
        } else {
            linf.AddMember("estimated", "false", allocator);
        }
        linf.AddMember("min", this->l_inf.min, allocator);
        linf.AddMember("min", this->l_inf.max, allocator);
        linf.AddMember("phase", this->l_inf.phase, allocator);
        parameters.AddMember("linf", linf, allocator);

        growth.AddMember("parameters", parameters, allocator);

        //empirical weight at age
        if (this->has_emprical_weight) {
            rapidjson::Value empirical_weight_at_age(rapidjson::kArrayType);

            if (this->weight_at_catch_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_catch_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_catch_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "catch_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_females, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_survey_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_survey_males, 3, nyears, nseasons,
                        nages);
                ewa.AddMember("data_object_type",
                        "survey_empirical_weight_at_age", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_season_start_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_season_start_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_season_start", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_females.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_females, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "females", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }

            if (this->weight_at_spawning_males.size() != 0) {
                rapidjson::Value ewa(rapidjson::kObjectType);
                rapidjson::Value values(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, values,
                        this->weight_at_spawning_males, 2, nyears, nages,
                        nages);
                ewa.AddMember("data_object_type",
                        "empirical_weight_at_age_spawning", allocator);
                ewa.AddMember("sex", "males", allocator);
                ewa.AddMember("units", "KG", allocator);
                ewa.AddMember("missing_values", "-999", allocator);
                ewa.AddMember("values", values, allocator);
                empirical_weight_at_age.PushBack(ewa, allocator);
            }
            growth.AddMember("empirical_weight_at_age", empirical_weight_at_age,
                    allocator);

        }
        grth.PushBack(growth, allocator);

    }

    static std::map<int, VonBertalanffyModified*> initialized_models;
    typedef typename std::map<int, VonBertalanffyModified*>::iterator model_iterator;
};

std::map<int, VonBertalanffyModified*> VonBertalanffyModified::initialized_models;

class Area : public MASSubModel {
protected:

public:

    static int id_g;
    int id;
    std::string name;

    Area() {

        id = Area::id_g++;
        Area::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~Area() {
    }

    virtual void AddToMAS(mas::Information<double> &info) {

        atl::intrusive_ptr<mas::Area<double> > area = new mas::Area<double>();
        mas::Area<double> *a = area.get();
        a->id = this->id;
        info.areas[a->id] = area;
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value area(rapidjson::kObjectType);

        area.AddMember("id", this->id, allocator);
        document.AddMember("area", area, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &a, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value area(rapidjson::kObjectType);

        area.AddMember("id", this->id, allocator);
        a.PushBack(area, allocator);
    }

    static std::map<int, Area*> initialized_models;
    typedef typename std::map<int, Area*>::iterator model_iterator;
};

std::map<int, Area*> Area::initialized_models;

int Area::id_g = 1;

/**
 * Movement
 */

class Movement : public MASSubModel {
protected:

public:
    static int id_g;
    int id;
    Rcpp::NumericVector connectivity_males;
    Rcpp::NumericVector connectivity_females;
    Rcpp::NumericVector connectivity_recruits;

    bool estimate_movement = false;
    int movement_phase = 1;

    Movement() {

        this->id = Movement::id_g++;
        Movement::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~Movement() {
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        if (this->connectivity_females.size() != std::pow(Area::id_g - 1, 2)
                || this->connectivity_males.size()
                != std::pow(Area::id_g - 1, 2)
                || this->connectivity_recruits.size()
                != std::pow(Area::id_g - 1, 2)) {
            std::cout << "MAS Error: Movement connectivity not equal to "
                    << std::pow(Area::id_g - 1, 2) << "\n";
            info.valid_configuration = false;
        } else {
            typedef typename mas::VariableTrait<double>::variable variable;

            atl::intrusive_ptr<mas::Movement<double> > movement =
                    new mas::Movement<double>();
            mas::Movement<double> *m = movement.get();
            m->id = this->id;
            int k = 0;
            m->male_connectivity.resize(info.nseasons);
            m->female_connectivity.resize(info.nseasons);
            m->recruit_connectivity.resize(info.nseasons);

            for (int s = 0; s < info.nseasons; s++) {
                m->male_connectivity[s].resize(Area::id_g - 1);
                m->female_connectivity[s].resize(Area::id_g - 1);
                m->recruit_connectivity[s].resize(Area::id_g - 1);

                for (int i = 0; i < Area::id_g - 1; i++) {
                    m->male_connectivity[s][i].resize(Area::id_g - 1);
                    m->female_connectivity[s][i].resize(Area::id_g - 1);
                    m->recruit_connectivity[s][i].resize(Area::id_g - 1);
                    for (int j = 0; j < Area::id_g - 1; j++) {
                        k = (i * (Area::id_g - 1)) + j;
                        m->male_connectivity[s][i][j] = variable(
                                this->connectivity_males[k]);
                        m->female_connectivity[s][i][j] = variable(
                                this->connectivity_females[k]);
                        m->recruit_connectivity[s][i][j] = variable(
                                this->connectivity_recruits[k]);
                        if (this->estimate_movement) {
                            m->Register(m->male_connectivity[s][i][j]);
                            m->Register(m->female_connectivity[s][i][j]);
                            m->Register(m->recruit_connectivity[s][i][j]);

                        }
                    }

                }
            }
            info.movement_models[m->id] = movement;
        }
    }

    void ExtractFromMAS(mas::Information<double> &info) {
        typename mas::Information<double>::movement_model_iterator mit;
        mit = info.movement_models.find(this->id);

        if (mit != info.movement_models.end()) {
            atl::intrusive_ptr<mas::Movement<double> > movement = (*mit).second;
            mas::Movement<double> *m = movement.get();
            int k = 0;
            for (int s = 0; s < info.nseasons; s++) {
                for (int i = 0; i < Area::id_g; i++) {
                    for (int j = 0; j < Area::id_g; j++) {

                        this->connectivity_males[k] =
                                m->male_connectivity[s][i][j].GetValue();
                        this->connectivity_females[k] =
                                m->female_connectivity[s][i][j].GetValue();
                        this->connectivity_recruits[k++] =
                                m->recruit_connectivity[s][i][j].GetValue();
                    }
                }
            }
        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value movement(rapidjson::kObjectType);

        movement.AddMember("id", this->id, allocator);
        if (this->estimate_movement) {
            movement.AddMember("estimated", "true", allocator);
        } else {
            movement.AddMember("estimated", "false", allocator);
        }

        rapidjson::Value recruits(rapidjson::kArrayType);
        rapidjson::Value females(rapidjson::kArrayType);
        rapidjson::Value males(rapidjson::kArrayType);
        if (this->connectivity_recruits.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for recruits vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, recruits,
                this->connectivity_recruits, 3, 1, nseasons, nareas);
        movement.AddMember("recruits", recruits, allocator);
        if (this->connectivity_females.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for females vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, females,
                this->connectivity_females, 3, 1, nseasons, nareas);
        movement.AddMember("female", females, allocator);

        if (this->connectivity_males.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for males vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, males,
                this->connectivity_males, 3, 1, nseasons, nareas);
        movement.AddMember("male", males, allocator);
        document.AddMember("movement", movement, allocator);

    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &move, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value movement(rapidjson::kObjectType);

        movement.AddMember("id", this->id, allocator);
        if (this->estimate_movement) {
            movement.AddMember("estimated", "true", allocator);
        } else {
            movement.AddMember("estimated", "false", allocator);
        }

        rapidjson::Value recruits(rapidjson::kArrayType);
        rapidjson::Value females(rapidjson::kArrayType);
        rapidjson::Value males(rapidjson::kArrayType);
        if (this->connectivity_recruits.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for recruits vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, recruits,
                this->connectivity_recruits, 3, 1, nseasons, nareas);
        movement.AddMember("recruits", recruits, allocator);
        if (this->connectivity_females.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for females vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, females,
                this->connectivity_females, 3, 1, nseasons, nareas);
        movement.AddMember("female", females, allocator);

        if (this->connectivity_males.size() != nseasons * nareas) {
            std::cout
                    << "MAS Warning: Movement connectivity matrix for males vector not equal to "
                    << nseasons * nareas << "." << std::endl;
        }
        MASSubModel::GenerateArrayObject(document, males,
                this->connectivity_males, 3, 1, nseasons, nareas);
        movement.AddMember("male", males, allocator);
        move.PushBack(movement, allocator);
    }

    static std::map<int, Movement*> initialized_models;
    typedef typename std::map<int, Movement*>::iterator model_iterator;
};

std::map<int, Movement*> Movement::initialized_models;

int Movement::id_g = 1;

class Maturity {
public:
    static int id_g;
    int id;
    Rcpp::NumericVector values;

    Maturity() {

        this->id = Maturity::id_g++;
        Maturity::initialized_models[this->id] = this;
    }

    virtual ~Maturity() {
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &mat, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value maturity(rapidjson::kObjectType);

        maturity.AddMember("id", this->id, allocator);
        rapidjson::Value values(rapidjson::kArrayType);

        for (int i = 0; i < this->values.size(); i++) {
            values.PushBack(this->values[i], allocator);
        }

        maturity.AddMember("values", values, allocator);
        mat.PushBack(maturity, allocator);
    }

    static std::map<int, Maturity*> initialized_models;
    typedef typename std::map<int, Maturity*>::iterator model_iterator;
};

std::map<int, Maturity*> Maturity::initialized_models;
int Maturity::id_g = 1;

class Population : public MASSubModel {
protected:
public:

    static int id_g;
    std::vector<std::pair<int, int> > movement; //id, year, connectivity
    std::vector<std::pair<int, int> > natural_mortality_males; //area,id
    std::vector<std::pair<int, int> > natural_mortality_females; //area,id
    std::map<int, int> intial_deviations_males; //area,id
    std::map<int, int> intial_deviations_females; //area,id
    std::vector<triple<int, int, int> > recruitment; //area,id
    int growth = -999;
    typedef typename std::map<int, int>::iterator deviations_iterator;
public:

    int id;
    double sex_ratio = 0.5;
    double spawning_season_offset = 0.35;

    //not exposed
    std::map<int, int> maturity_males;
    std::map<int, int> maturity_females;
    typedef typename std::map<int, int>::iterator maturity_iterator;

    Population() {

        this->id = Population::id_g++;
        Population::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~Population() {
    }

    void SetGrowth(int id) {

        this->growth = id;
        std::cout << "Setting growth " << this->growth << "\n";
    }

    void AddMovement(int id, int year) {

        this->movement.push_back(std::make_pair(id, year));
    }

    void AddMaturity(int id, int area, std::string sex) {

        std::locale loc;
        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        switch (GetSexType(sex)) {

            case MALE:
                this->maturity_males[area] = id;
                break;
            case FEMALE:
                this->maturity_females[area] = id;
                break;

            case UNDIFFERENTIATED:
                this->maturity_males[area] = id;
                this->maturity_females[area] = id;

                break;
            default:
                std::cout << "MAS Error: unknown sex type \"" << sex
                        << "\" for maturity input.\n";
        }
    }

    void AddNaturalMortality(int id, int area, const std::string &sex) {
        std::locale loc;
        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        switch (GetSexType(sex)) {

            case MALE:
                this->natural_mortality_males.push_back(std::make_pair(id, area));
                break;
            case FEMALE:
                this->natural_mortality_females.push_back(std::make_pair(id, area));
                break;

            case UNDIFFERENTIATED:
                this->natural_mortality_males.push_back(std::make_pair(id, area));
                this->natural_mortality_females.push_back(std::make_pair(id, area));

                break;
            default:
                std::cout << "MAS Error: unknown sex type \"" << sex
                        << "\" for natural mortality input.\n";
        }
    }

    void AddRecruitment(int id, int season, int area) {
        this->recruitment.push_back(triple<int, int, int>(id, season, area));
    }

    void SetInitialDeviations(int devs, int area, std::string sex) {

        std::locale loc;
        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        switch (GetSexType(sex)) {

            case MALE:
                this->intial_deviations_males[area] = devs;
                break;
            case FEMALE:
                this->intial_deviations_females[area] = devs;
                break;

            case UNDIFFERENTIATED:
                this->intial_deviations_males[area] = devs;
                this->intial_deviations_females[area] = devs;

                break;
            default:
                std::cout << "MAS Error: unknown sex type \"" << sex
                        << "\" for maturity input.\n";
        }
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        typedef typename mas::VariableTrait<double>::variable variable;
        atl::intrusive_ptr<mas::Population<double> > population =
                new mas::Population<double>();
        mas::Population<double> *pop = population.get();

        pop->id = this->id;
        pop->female_fraction_value = this->sex_ratio;
        pop->growth_id = this->growth;
        maturity_iterator mit;
        /**
         * Maturity
         */
        for (mit = this->maturity_females.begin();
                mit != this->maturity_females.end(); ++mit) {
            typename Maturity::model_iterator it;
            it = Maturity::initialized_models.find((*mit).first);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &v = (*it).second->values;
                std::vector<double> values(v.size());
                for (int i = 0; i < v.size(); i++) {
                    values[i] = v[i];
                }
                pop->maturity_models[(*mit).second][mas::FEMALE] = values;
            }
        }

        for (mit = this->maturity_males.begin();
                mit != this->maturity_males.end(); ++mit) {
            typename Maturity::model_iterator it;
            it = Maturity::initialized_models.find((*mit).first);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &v = (*it).second->values;
                std::vector<double> values(v.size());
                for (int i = 0; i < v.size(); i++) {
                    values[i] = v[i];
                }
                pop->maturity_models[(*mit).second][mas::MALE] = values;
            }
        }

        for (int i = 0; i < this->movement.size(); i++) {
            pop->movement_models_ids[movement[i].second] = movement[i].first;
        }

        /*
         * Initial deviations
         */
        deviations_iterator dev_it;
        for (dev_it = this->intial_deviations_females.begin();
                dev_it != this->intial_deviations_females.end(); ++dev_it) {
            typename InitialDeviations::model_iterator dit;

            dit = InitialDeviations::initialized_models.find((*dev_it).second);
            if (dit != InitialDeviations::initialized_models.end()) {
                InitialDeviations *devs = (*dit).second;
                if (devs->values.size() != info.ages.size()) {
                    std::cout
                            << "MAS Error: Initial deviations vector size not equal to "
                            << info.ages.size() << "\n";
                    info.valid_configuration = false;
                } else {

                    std::vector<variable> &v =
                            pop->initial_deviations_females[(*dev_it).first].second;
                    v.resize(devs->values.size());
                    for (int i = 0; i < v.size(); i++) {
                        v[i] = variable(devs->values[i]);
                    }

                    if (devs->estimate) {
                        for (int i = 0; i < v.size(); i++) {
                            std::stringstream ss;
                            ss << "f_initial_deviations[" << i << "]_"
                                    << pop->id;
                            mas::VariableTrait<double>::SetName(v[i], ss.str());
                            pop->Register(v[i], devs->phase);
                        }
                    }

                }

            }
        }

        for (dev_it = this->intial_deviations_males.begin();
                dev_it != this->intial_deviations_males.end(); ++dev_it) {
            typename InitialDeviations::model_iterator dit;

            dit = InitialDeviations::initialized_models.find((*dev_it).second);
            if (dit != InitialDeviations::initialized_models.end()) {
                InitialDeviations *devs = (*dit).second;
                if (devs->values.size() != info.ages.size()) {
                    std::cout
                            << "MAS Error: Initial deviations vector size not equal to "
                            << info.ages.size() << "\n";
                    info.valid_configuration = false;
                } else {

                    std::vector<variable> &v =
                            pop->initial_deviations_males[(*dev_it).first].second;
                    v.resize(devs->values.size());
                    for (int i = 0; i < v.size(); i++) {
                        v[i] = variable(devs->values[i]);
                    }

                    if (devs->estimate) {
                        for (int i = 0; i < v.size(); i++) {
                            std::stringstream ss;
                            ss << "m_initial_deviations[" << i << "]_"
                                    << pop->id;
                            mas::VariableTrait<double>::SetName(v[i], ss.str());
                            pop->Register(v[i], devs->phase);
                        }
                    }

                }

            }
        }

        /*
         * Recruitment
         */
        for (int i = 0; i < this->recruitment.size(); i++) {
            pop->recruitment_ids[recruitment[i].third] =
                    this->recruitment[i].first;
        }

        std::cout << "growth model for this population is " << this->growth << "\n";
        /*
         * Growth
         */
        if (this->growth != -999) {
            pop->growth_id = growth;
        } else {
            std::cout << "growth id " << growth << "\n\n";
            std::cout << "MAS Error: No growth model set for Population "
                    << this->id << "\n";
            info.valid_configuration = false;
        }

        /*
         * Natural Mortality
         */

        for (int i = 0; i < this->natural_mortality_males.size(); i++) {
            pop->male_natural_mortality_ids[this->natural_mortality_males[i].second] =
                    this->natural_mortality_males[i].first;
        }
        for (int i = 0; i < this->natural_mortality_females.size(); i++) {

            pop->female_natural_mortality_ids[this->natural_mortality_females[i].second] =
                    this->natural_mortality_females[i].first;
        }

        info.populations[this->id] = population;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value population(rapidjson::kObjectType);
        population.AddMember("id", this->id, allocator);
        population.AddMember("sex_ratio", this->sex_ratio, allocator);

        rapidjson::Value movement(rapidjson::kArrayType);
        for (int i = 0; i < this->movement.size(); i++) {

            rapidjson::Value movement_entry(rapidjson::kObjectType);
            movement_entry.AddMember("year", this->movement[i].second,
                    allocator);
            movement_entry.AddMember("id", this->movement[i].first, allocator);
            movement.PushBack(movement_entry, allocator);
        }

        population.AddMember("movement", movement, allocator);

        rapidjson::Value maturity(rapidjson::kArrayType);
        //females
        maturity_iterator mit;
        for (mit = this->maturity_females.begin();
                mit != this->maturity_females.end(); ++mit) {
            rapidjson::Value maturity_entry(rapidjson::kObjectType);
            maturity_entry.AddMember("sex", "females", allocator);
            maturity_entry.AddMember("area", (*mit).first, allocator);

            Maturity::model_iterator it = Maturity::initialized_models.find(
                    (*mit).second);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &values = (*it).second->values;
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < values.size(); i++) {
                    vals.PushBack(values[i], allocator);
                }
                maturity_entry.AddMember("values", vals, allocator);
            }
            maturity.PushBack(maturity_entry, allocator);
        }

        for (mit = this->maturity_males.begin();
                mit != this->maturity_males.end(); ++mit) {
            rapidjson::Value maturity_entry(rapidjson::kObjectType);
            maturity_entry.AddMember("sex", "females", allocator);
            maturity_entry.AddMember("area", (*mit).first, allocator);

            Maturity::model_iterator it = Maturity::initialized_models.find(
                    (*mit).second);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &values = (*it).second->values;
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < values.size(); i++) {
                    vals.PushBack(values[i], allocator);
                }
                maturity_entry.AddMember("values", vals, allocator);
            }
            maturity.PushBack(maturity_entry, allocator);
        }
        population.AddMember("maturity", maturity, allocator);

        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value natural_mortality(rapidjson::kArrayType);
        //females
        for (int i = 0; i < this->natural_mortality_females.size(); i++) {
            rapidjson::Value natural_mortality_entry(rapidjson::kObjectType);
            natural_mortality_entry.AddMember("sex", "female", allocator);
            natural_mortality_entry.AddMember("id",
                    natural_mortality_females[i].first, allocator);
            natural_mortality_entry.AddMember("area",
                    natural_mortality_females[i].second, allocator);
            natural_mortality.PushBack(natural_mortality_entry, allocator);
        }
        //males
        for (int i = 0; i < this->natural_mortality_males.size(); i++) {
            rapidjson::Value natural_mortality_entry(rapidjson::kObjectType);
            natural_mortality_entry.AddMember("sex", "male", allocator);
            natural_mortality_entry.AddMember("id",
                    natural_mortality_males[i].first, allocator);
            natural_mortality_entry.AddMember("area",
                    natural_mortality_males[i].second, allocator);
            natural_mortality.PushBack(natural_mortality_entry, allocator);
        }

        parameters.AddMember("natural_mortality", natural_mortality, allocator);

        rapidjson::Value recruitment(rapidjson::kArrayType);
        for (int i = 0; i < this->recruitment.size(); i++) {
            rapidjson::Value recruitment_entry(rapidjson::kObjectType);
            recruitment_entry.AddMember("id", this->recruitment[i].first,
                    allocator);
            recruitment_entry.AddMember("season", this->recruitment[i].second,
                    allocator);
            recruitment_entry.AddMember("area", this->recruitment[i].third,
                    allocator);
            recruitment.PushBack(recruitment_entry, allocator);
        }
        parameters.AddMember("recruitment", recruitment, allocator);

        rapidjson::Value growth(rapidjson::kObjectType);
        growth.AddMember("id", this->growth, allocator);
        parameters.AddMember("growth", growth, allocator);

        rapidjson::Value initial_deviations(rapidjson::kArrayType);
        deviations_iterator idit;
        for (idit = this->intial_deviations_females.begin();
                idit != this->intial_deviations_females.end(); ++idit) {
            rapidjson::Value initial_devs_female(rapidjson::kObjectType);
            initial_devs_female.AddMember("sex", "female", allocator);
            initial_devs_female.AddMember("area", (*idit).first, allocator);

            InitialDeviations::model_iterator it =
                    InitialDeviations::initialized_models.find((*idit).second);
            if (it != InitialDeviations::initialized_models.end()) {
                InitialDeviations *idevs = it->second;
                if (idevs->estimate == true) {
                    initial_devs_female.AddMember("estimated", "true",
                            allocator);
                } else {
                    initial_devs_female.AddMember("estimated", "false",
                            allocator);
                }
                initial_devs_female.AddMember("random_effect", "false",
                        allocator);
                initial_devs_female.AddMember("id", idevs->id, allocator);
                initial_devs_female.AddMember("phase", idevs->phase, allocator);
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < idevs->values.size(); i++) {
                    vals.PushBack(idevs->values[i], allocator);
                }
                initial_devs_female.AddMember("values", vals, allocator);
            }
            initial_deviations.PushBack(initial_devs_female, allocator);
        }

        for (idit = this->intial_deviations_males.begin();
                idit != this->intial_deviations_males.end(); ++idit) {
            rapidjson::Value initial_devs_male(rapidjson::kObjectType);
            initial_devs_male.AddMember("sex", "male", allocator);
            initial_devs_male.AddMember("area", (*idit).first, allocator);

            InitialDeviations::model_iterator it =
                    InitialDeviations::initialized_models.find((*idit).second);
            if (it != InitialDeviations::initialized_models.end()) {
                InitialDeviations *idevs = it->second;
                if (idevs->estimate == true) {
                    initial_devs_male.AddMember("estimated", "true", allocator);
                } else {
                    initial_devs_male.AddMember("estimated", "false",
                            allocator);
                }
                initial_devs_male.AddMember("random_effect", "false",
                        allocator);
                initial_devs_male.AddMember("id", idevs->id, allocator);
                initial_devs_male.AddMember("phase", idevs->phase, allocator);
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < idevs->values.size(); i++) {
                    vals.PushBack(idevs->values[i], allocator);
                }
                initial_devs_male.AddMember("values", vals, allocator);
            }
            initial_deviations.PushBack(initial_devs_male, allocator);
        }
        parameters.AddMember("initial_deviations", initial_deviations,
                allocator);

        population.AddMember("parameters", parameters, allocator);
        document.AddMember("population", population, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &pops, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value population(rapidjson::kObjectType);
        population.AddMember("id", this->id, allocator);
        population.AddMember("sex_ratio", this->sex_ratio, allocator);

        rapidjson::Value movement(rapidjson::kArrayType);
        for (int i = 0; i < this->movement.size(); i++) {

            rapidjson::Value movement_entry(rapidjson::kObjectType);
            movement_entry.AddMember("year", this->movement[i].second,
                    allocator);
            movement_entry.AddMember("id", this->movement[i].first, allocator);
            movement.PushBack(movement_entry, allocator);
        }

        population.AddMember("movement", movement, allocator);

        rapidjson::Value maturity(rapidjson::kArrayType);
        //females
        maturity_iterator mit;
        for (mit = this->maturity_females.begin();
                mit != this->maturity_females.end(); ++mit) {
            rapidjson::Value maturity_entry(rapidjson::kObjectType);
            maturity_entry.AddMember("sex", "females", allocator);
            maturity_entry.AddMember("area", (*mit).first, allocator);

            Maturity::model_iterator it = Maturity::initialized_models.find(
                    (*mit).second);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &values = (*it).second->values;
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < values.size(); i++) {
                    vals.PushBack(values[i], allocator);
                }
                maturity_entry.AddMember("values", vals, allocator);
            }
            maturity.PushBack(maturity_entry, allocator);
        }

        for (mit = this->maturity_males.begin();
                mit != this->maturity_males.end(); ++mit) {
            rapidjson::Value maturity_entry(rapidjson::kObjectType);
            maturity_entry.AddMember("sex", "females", allocator);
            maturity_entry.AddMember("area", (*mit).first, allocator);

            Maturity::model_iterator it = Maturity::initialized_models.find(
                    (*mit).second);
            if (it != Maturity::initialized_models.end()) {
                Rcpp::NumericVector &values = (*it).second->values;
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < values.size(); i++) {
                    vals.PushBack(values[i], allocator);
                }
                maturity_entry.AddMember("values", vals, allocator);
            }
            maturity.PushBack(maturity_entry, allocator);
        }
        population.AddMember("maturity", maturity, allocator);

        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value natural_mortality(rapidjson::kArrayType);
        //females
        for (int i = 0; i < this->natural_mortality_females.size(); i++) {
            rapidjson::Value natural_mortality_entry(rapidjson::kObjectType);
            natural_mortality_entry.AddMember("sex", "female", allocator);
            natural_mortality_entry.AddMember("id",
                    natural_mortality_females[i].first, allocator);
            natural_mortality_entry.AddMember("area",
                    natural_mortality_females[i].second, allocator);
            natural_mortality.PushBack(natural_mortality_entry, allocator);
        }
        //males
        for (int i = 0; i < this->natural_mortality_males.size(); i++) {
            rapidjson::Value natural_mortality_entry(rapidjson::kObjectType);
            natural_mortality_entry.AddMember("sex", "male", allocator);
            natural_mortality_entry.AddMember("id",
                    natural_mortality_males[i].first, allocator);
            natural_mortality_entry.AddMember("area",
                    natural_mortality_males[i].second, allocator);
            natural_mortality.PushBack(natural_mortality_entry, allocator);
        }

        parameters.AddMember("natural_mortality", natural_mortality, allocator);

        rapidjson::Value recruitment(rapidjson::kArrayType);
        for (int i = 0; i < this->recruitment.size(); i++) {
            rapidjson::Value recruitment_entry(rapidjson::kObjectType);
            recruitment_entry.AddMember("id", this->recruitment[i].first,
                    allocator);
            recruitment_entry.AddMember("season", this->recruitment[i].second,
                    allocator);
            recruitment_entry.AddMember("area", this->recruitment[i].third,
                    allocator);
            recruitment.PushBack(recruitment_entry, allocator);
        }
        parameters.AddMember("recruitment", recruitment, allocator);

        rapidjson::Value growth(rapidjson::kObjectType);
        growth.AddMember("id", this->growth, allocator);
        parameters.AddMember("growth", growth, allocator);

        rapidjson::Value initial_deviations(rapidjson::kArrayType);
        deviations_iterator idit;
        for (idit = this->intial_deviations_females.begin();
                idit != this->intial_deviations_females.end(); ++idit) {
            rapidjson::Value initial_devs_female(rapidjson::kObjectType);
            initial_devs_female.AddMember("sex", "female", allocator);
            initial_devs_female.AddMember("area", (*idit).first, allocator);

            InitialDeviations::model_iterator it =
                    InitialDeviations::initialized_models.find((*idit).second);
            if (it != InitialDeviations::initialized_models.end()) {
                InitialDeviations *idevs = it->second;
                if (idevs->estimate == true) {
                    initial_devs_female.AddMember("estimated", "true",
                            allocator);
                } else {
                    initial_devs_female.AddMember("estimated", "false",
                            allocator);
                }
                initial_devs_female.AddMember("random_effect", "false",
                        allocator);
                initial_devs_female.AddMember("id", idevs->id, allocator);
                initial_devs_female.AddMember("phase", idevs->phase, allocator);
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < idevs->values.size(); i++) {
                    vals.PushBack(idevs->values[i], allocator);
                }
                initial_devs_female.AddMember("values", vals, allocator);
            }
            initial_deviations.PushBack(initial_devs_female, allocator);
        }

        for (idit = this->intial_deviations_males.begin();
                idit != this->intial_deviations_males.end(); ++idit) {
            rapidjson::Value initial_devs_male(rapidjson::kObjectType);
            initial_devs_male.AddMember("sex", "male", allocator);
            initial_devs_male.AddMember("area", (*idit).first, allocator);

            InitialDeviations::model_iterator it =
                    InitialDeviations::initialized_models.find((*idit).second);
            if (it != InitialDeviations::initialized_models.end()) {
                InitialDeviations *idevs = it->second;
                if (idevs->estimate == true) {
                    initial_devs_male.AddMember("estimated", "true", allocator);
                } else {
                    initial_devs_male.AddMember("estimated", "false",
                            allocator);
                }
                initial_devs_male.AddMember("random_effect", "false",
                        allocator);
                initial_devs_male.AddMember("id", idevs->id, allocator);
                initial_devs_male.AddMember("phase", idevs->phase, allocator);
                rapidjson::Value vals(rapidjson::kArrayType);
                for (int i = 0; i < idevs->values.size(); i++) {
                    vals.PushBack(idevs->values[i], allocator);
                }
                initial_devs_male.AddMember("values", vals, allocator);
            }
            initial_deviations.PushBack(initial_devs_male, allocator);
        }
        parameters.AddMember("initial_deviations", initial_deviations,
                allocator);

        population.AddMember("parameters", parameters, allocator);
        pops.PushBack(population, allocator);
    }

    static std::map<int, Population*> initialized_models;
    typedef typename std::map<int, Population*>::iterator model_iterator;
};

std::map<int, Population*> Population::initialized_models;
int Population::id_g = 1;

/**
 * NLL Components
 */

class NLLBase : public MASSubModel {
public:

    static int id_g;
    int id_;
    static std::vector<NLLBase*> nll_submodels;
    typedef typename std::vector<NLLBase*>::iterator nll_iterator;

    virtual ~NLLBase() {

    }
};
int NLLBase::id_g = 1;
std::vector<NLLBase*> NLLBase::nll_submodels;

class Lognormal : public NLLBase {
    Rcpp::NumericVector lambdas;
    Rcpp::IntegerVector lambda_dimensions;
    bool has_lambdas = false;
public:
    int id;
    bool use_bias_correction = true;

    Lognormal() {

        this->id = NLLBase::id_g++;
        this->id_ = this->id;
        Lognormal::initialized_models[this->id] = this;
        //MASSubModel::submodels.push_back(this);
        NLLBase::nll_submodels.push_back(this);
    }

    virtual ~Lognormal() {
    }

    void SetLambdaValues(Rcpp::NumericVector lambdas,
            Rcpp::IntegerVector lambda_dimensions) {

        this->lambdas = lambdas;
        this->lambda_dimensions = lambda_dimensions;
        this->has_lambdas = true;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::DataObject<double> > lambda_data(
                new mas::DataObject<double>());
        atl::intrusive_ptr<mas::Lognormal<double> > ln = new mas::Lognormal<
                double>();
        mas::Lognormal<double> *nll = ln.get();
        ln->use_bias_correction = this->use_bias_correction;
        nll->lambda = lambda_data;
        nll->id = this->id;
        if (this->has_lambdas) {
            nll->lambda = new mas::DataObject<double>();
            mas::DataObject<double> *d = nll->lambda.get();
            int prod = 1.0;
            for (int i = 0; i < this->lambda_dimensions.size(); i++) {
                prod *= this->lambdas[i];
            }

            if (this->lambdas.size() != prod) {
                std::cout
                        << "MAS Error: Lognormal lambda vector not equal to dimensions product\n";
                info.valid_configuration = false;
            }

            switch (lambdas.size()) {
                case 1:
                    d->imax = this->lambda_dimensions[0];
                    break;
                case 2:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    break;
                case 3:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    d->kmax = this->lambda_dimensions[2];
                    break;
            }
            for (int i = 0; i < this->lambdas.size(); i++) {

                d->data.push_back(this->lambdas[i]);
            }
        }

        info.likelihood_components[nll->id] = ln;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "lognormal", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 2) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        2, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"lognormal\" \"lambdas\" entry.\n";
            }
        }
        document.AddMember("likelihood_component", likelihood, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nll, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "lognormal", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 2) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        2, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"lognormal\" \"lambdas\" entry.\n";
            }
        }
        nll.PushBack(likelihood, allocator);
    }

    static std::map<int, Lognormal*> initialized_models;
    typedef typename std::map<int, Lognormal*>::iterator model_iterator;
};

std::map<int, Lognormal*> Lognormal::initialized_models;

class DirichletMultinomial : public NLLBase {
    Rcpp::NumericVector lambdas;
    Rcpp::IntegerVector lambda_dimensions;
    bool has_lambdas = false;
public:
    int id;
    Parameter beta;

    DirichletMultinomial() {

        this->id = NLLBase::id_g++;
        this->id_ = this->id;

        DirichletMultinomial::initialized_models[this->id] = this;
        //MASSubModel::submodels.push_back(this);
        NLLBase::nll_submodels.push_back(this);
    }

    virtual ~DirichletMultinomial() {
    }

    void SetLambdaValues(Rcpp::NumericVector lambdas,
            Rcpp::IntegerVector lambda_dimensions) {

        this->lambdas = lambdas;
        this->lambda_dimensions = lambda_dimensions;
        this->has_lambdas = true;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::DataObject<double> > lambda_data(
                new mas::DataObject<double>());
        atl::intrusive_ptr<mas::DirichletMultinomial<double> > ln =
                new mas::DirichletMultinomial<double>();
        mas::DirichletMultinomial<double> *nll = ln.get();
        nll->lambda = lambda_data;
        nll->id = this->id;
        if (this->has_lambdas) {
            nll->lambda = new mas::DataObject<double>();
            mas::DataObject<double> *d = nll->lambda.get();
            int prod = 1.0;
            for (int i = 0; i < this->lambda_dimensions.size(); i++) {
                prod *= this->lambdas[i];
            }

            if (this->lambdas.size() != prod) {
                std::cout
                        << "MAS Error: DirichletMultinomial lambda vector not equal to dimensions product\n";
                info.valid_configuration = false;
            }

            switch (lambdas.size()) {
                case 1:
                    d->imax = this->lambda_dimensions[0];
                    break;
                case 2:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    break;
                case 3:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    d->kmax = this->lambda_dimensions[2];
                    break;
            }
            for (int i = 0; i < this->lambdas.size(); i++) {

                d->data.push_back(this->lambdas[i]);
            }
        }

        info.likelihood_components[nll->id] = ln;
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "dirichlet_multinomial", allocator);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
        } else {
            beta.AddMember("estimated", "false", allocator);
        }
        beta.AddMember("min", this->beta.min, allocator);
        beta.AddMember("max", this->beta.max, allocator);
        beta.AddMember("phase", this->beta.phase, allocator);
        parameters.AddMember("beta", beta, allocator);
        likelihood.AddMember("parameters", parameters, allocator);
        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"dirichlet_multinomial\" \"lambdas\" entry.\n";
            }
        }
        document.AddMember("likelihood_component", likelihood, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nll, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "dirichlet_multinomial", allocator);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
        } else {
            beta.AddMember("estimated", "false", allocator);
        }
        beta.AddMember("min", this->beta.min, allocator);
        beta.AddMember("max", this->beta.max, allocator);
        beta.AddMember("phase", this->beta.phase, allocator);
        parameters.AddMember("beta", beta, allocator);
        likelihood.AddMember("parameters", parameters, allocator);
        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"dirichlet_multinomial\" \"lambdas\" entry.\n";
            }
        }
        nll.PushBack(likelihood, allocator);
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    static std::map<int, DirichletMultinomial*> initialized_models;
    typedef typename std::map<int, DirichletMultinomial*>::iterator model_iterator;
};

std::map<int, DirichletMultinomial*> DirichletMultinomial::initialized_models;

class DirichletMultinomialRobust : public NLLBase {
    Rcpp::NumericVector lambdas;
    Rcpp::IntegerVector lambda_dimensions;
    bool has_lambdas = false;
public:
    int id;
    Parameter beta;

    DirichletMultinomialRobust() {

        this->id = NLLBase::id_g++;
        this->id_ = this->id;

        DirichletMultinomialRobust::initialized_models[this->id] = this;
        //MASSubModel::submodels.push_back(this);
        NLLBase::nll_submodels.push_back(this);
    }

    virtual ~DirichletMultinomialRobust() {
    }

    void SetLambdaValues(Rcpp::NumericVector lambdas,
            Rcpp::IntegerVector lambda_dimensions) {

        this->lambdas = lambdas;
        this->lambda_dimensions = lambda_dimensions;
        this->has_lambdas = true;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::DataObject<double> > lambda_data(
                new mas::DataObject<double>());

        atl::intrusive_ptr<mas::DirichletMultinomialRobust<double> > ln =
                new mas::DirichletMultinomialRobust<double>();
        mas::DirichletMultinomialRobust<double> *nll = ln.get();
        nll->lambda = lambda_data;
        nll->id = this->id;
        if (this->has_lambdas) {
            nll->lambda = new mas::DataObject<double>();
            mas::DataObject<double> *d = nll->lambda.get();
            int prod = 1.0;
            for (int i = 0; i < this->lambda_dimensions.size(); i++) {
                prod *= this->lambdas[i];
            }

            if (this->lambdas.size() != prod) {
                std::cout
                        << "MAS Error: DirichletMultinomialRobust lambda vector not equal to dimensions product\n";
                info.valid_configuration = false;
            }

            switch (lambdas.size()) {
                case 1:
                    d->imax = this->lambda_dimensions[0];
                    break;
                case 2:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    break;
                case 3:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    d->kmax = this->lambda_dimensions[2];
                    break;
            }
            for (int i = 0; i < this->lambdas.size(); i++) {

                d->data.push_back(this->lambdas[i]);
            }
        }

        info.likelihood_components[nll->id] = ln;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "dirichlet_multinomial_robust",
                allocator);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
        } else {
            beta.AddMember("estimated", "false", allocator);
        }
        beta.AddMember("min", this->beta.min, allocator);
        beta.AddMember("max", this->beta.max, allocator);
        beta.AddMember("phase", this->beta.phase, allocator);
        parameters.AddMember("beta", beta, allocator);
        likelihood.AddMember("parameters", parameters, allocator);
        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"dirichlet_multinomial_robust\" \"lambdas\" entry.\n";
            }
        }
        document.AddMember("likelihood_component", likelihood, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nll, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "dirichlet_multinomial_robust",
                allocator);
        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value beta(rapidjson::kObjectType);
        beta.AddMember("value", this->beta.value, allocator);
        if (this->beta.estimated) {
            beta.AddMember("estimated", "true", allocator);
        } else {
            beta.AddMember("estimated", "false", allocator);
        }
        beta.AddMember("min", this->beta.min, allocator);
        beta.AddMember("max", this->beta.max, allocator);
        beta.AddMember("phase", this->beta.phase, allocator);
        parameters.AddMember("beta", beta, allocator);
        likelihood.AddMember("parameters", parameters, allocator);
        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"dirichlet_multinomial_robust\" \"lambdas\" entry.\n";
            }
        }
        nll.PushBack(likelihood, allocator);
    }

    static std::map<int, DirichletMultinomialRobust*> initialized_models;
    typedef typename std::map<int, DirichletMultinomialRobust*>::iterator model_iterator;
};

std::map<int, DirichletMultinomialRobust*> DirichletMultinomialRobust::initialized_models;

class Multinomial : public NLLBase {
    Rcpp::NumericVector lambdas;
    Rcpp::IntegerVector lambda_dimensions;
    bool has_lambdas = false;
public:
    int id;

    Multinomial() {

        this->id = NLLBase::id_g++;
        this->id_ = this->id;

        Multinomial::initialized_models[this->id] = this;
        //MASSubModel::submodels.push_back(this);
        NLLBase::nll_submodels.push_back(this);
    }

    virtual ~Multinomial() {
    }

    void SetLambdaValues(Rcpp::NumericVector lambdas,
            Rcpp::IntegerVector lambda_dimensions) {

        this->lambdas = lambdas;
        this->lambda_dimensions = lambda_dimensions;
        this->has_lambdas = true;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::Multinomial<double> > ln = new mas::Multinomial<
                double>();
        mas::Multinomial<double> *nll = ln.get();
        atl::intrusive_ptr<mas::DataObject<double> > lambda_data(
                new mas::DataObject<double>());
        nll->lambda = lambda_data;

        nll->id = this->id;
        if (this->has_lambdas) {
            nll->lambda = new mas::DataObject<double>();
            mas::DataObject<double> *d = nll->lambda.get();
            int prod = 1.0;
            for (int i = 0; i < this->lambda_dimensions.size(); i++) {
                prod *= this->lambdas[i];
            }

            if (this->lambdas.size() != prod) {
                std::cout
                        << "MAS Error: Multinomial lambda vector not equal to dimensions product\n";
                info.valid_configuration = false;
            }

            switch (lambdas.size()) {
                case 1:
                    d->imax = this->lambda_dimensions[0];
                    break;
                case 2:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    break;
                case 3:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    d->kmax = this->lambda_dimensions[2];
                    break;
            }
            for (int i = 0; i < this->lambdas.size(); i++) {

                d->data.push_back(this->lambdas[i]);
            }
        }

        info.likelihood_components[nll->id] = ln;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "multinomial", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"multinomial\" \"lambdas\" entry.\n";
            }
        }
        document.AddMember("likelihood_component", likelihood, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nll, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "multinomial", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"multinomial\" \"lambdas\" entry.\n";
            }
        }
        nll.PushBack(likelihood, allocator);
    }

    static std::map<int, Multinomial*> initialized_models;
    typedef typename std::map<int, Multinomial*>::iterator model_iterator;
};

std::map<int, Multinomial*> Multinomial::initialized_models;

class MultinomialRobust : public NLLBase {
    Rcpp::NumericVector lambdas;
    Rcpp::IntegerVector lambda_dimensions;
    bool has_lambdas = false;
public:
    int id;

    MultinomialRobust() {

        this->id = NLLBase::id_g++;
        this->id_ = this->id;

        MultinomialRobust::initialized_models[this->id] = this;
        //MASSubModel::submodels.push_back(this);
        NLLBase::nll_submodels.push_back(this);
    }

    virtual ~MultinomialRobust() {
    }

    void SetLambdaValues(Rcpp::NumericVector lambdas,
            Rcpp::IntegerVector lambda_dimensions) {

        this->lambdas = lambdas;
        this->lambda_dimensions = lambda_dimensions;
        this->has_lambdas = true;
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        atl::intrusive_ptr<mas::MultinomialRobust<double> > ln =
                new mas::MultinomialRobust<double>();
        mas::MultinomialRobust<double> *nll = ln.get();
        atl::intrusive_ptr<mas::DataObject<double> > lambda_data(
                new mas::DataObject<double>());
        nll->lambda = lambda_data;
        nll->id = this->id;
        if (this->has_lambdas) {
            nll->lambda = new mas::DataObject<double>();
            mas::DataObject<double> *d = nll->lambda.get();
            int prod = 1.0;
            for (int i = 0; i < this->lambda_dimensions.size(); i++) {
                prod *= this->lambdas[i];
            }

            if (this->lambdas.size() != prod) {
                std::cout
                        << "MAS Error: MultinomialRobust lambda vector not equal to dimensions product\n";
                info.valid_configuration = false;
            }

            switch (lambdas.size()) {
                case 1:
                    d->imax = this->lambda_dimensions[0];
                    break;
                case 2:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    break;
                case 3:
                    d->imax = this->lambda_dimensions[0];
                    d->jmax = this->lambda_dimensions[1];
                    d->kmax = this->lambda_dimensions[2];
                    break;
            }
            for (int i = 0; i < this->lambdas.size(); i++) {

                d->data.push_back(this->lambdas[i]);
            }
        }

        info.likelihood_components[nll->id] = ln;
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "multinomial_robust", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"multinomial_robust\" \"lambdas\" entry.\n";
            }
        }
        document.AddMember("likelihood_component", likelihood, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &nll, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value likelihood(rapidjson::kObjectType);
        likelihood.AddMember("id", this->id, allocator);
        likelihood.AddMember("model", "multinomial_robust", allocator);

        if (this->lambdas.size()) {
            if (this->lambda_dimensions.size() == 3) {
                rapidjson::Value lambdas(rapidjson::kObjectType);
                rapidjson::Value vals(rapidjson::kArrayType);
                MASSubModel::GenerateArrayObject(document, vals, this->lambdas,
                        3, this->lambda_dimensions[0],
                        this->lambda_dimensions[1], this->lambda_dimensions[2]);
                lambdas.AddMember("values", vals, allocator);
                likelihood.AddMember("lambdas", lambdas, allocator);
            } else {
                std::cout
                        << "MAS Warning: Incorrect dimensions for \"multinomial_robust\" \"lambdas\" entry.\n";
            }
        }
        nll.PushBack(likelihood, allocator);
    }

    static std::map<int, MultinomialRobust*> initialized_models;
    typedef typename std::map<int, MultinomialRobust*>::iterator model_iterator;
};

std::map<int, MultinomialRobust*> MultinomialRobust::initialized_models;

class IndexData : public MASSubModel {
public:
    static int id_g;
    bool is_abundance = false;

    IndexData(bool add_to_list = true) {

        if (add_to_list) {
            this->id = IndexData::id_g++;
            IndexData::initialized_models[this->id] = this;
        }
    }

    virtual ~IndexData() {
    }

    MASObjectType sub_model_type = DATAOBJECT;
    int id;
    Rcpp::NumericVector data;
    Rcpp::NumericVector error;
    std::string sex;
    double missing_values = -999;

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
    }

    static std::map<int, IndexData*> initialized_models;
    typedef typename std::map<int, IndexData*>::iterator model_iterator;
};

std::map<int, IndexData*> IndexData::initialized_models;
int IndexData::id_g = 1;

class AgeCompData : public MASSubModel {
public:
    static int id_g;

    AgeCompData(bool add_to_list = true) {

        if (add_to_list) {
            this->id = AgeCompData::id_g++;
            AgeCompData::initialized_models[this->id] = this;
        }
    }

    virtual ~AgeCompData() {
    }

    MASObjectType sub_model_type = DATAOBJECT;
    int id;
    Rcpp::NumericVector data;
    Rcpp::NumericVector sample_size;
    std::string sex;
    double missing_values = -999;

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
    }

    static std::map<int, AgeCompData*> initialized_models;
    typedef typename std::map<int, AgeCompData*>::iterator model_iterator;
};

std::map<int, AgeCompData*> AgeCompData::initialized_models;
int AgeCompData::id_g = 1;

class LengthCompData : public MASSubModel {
public:
    static int id_g;

    LengthCompData() {

        this->id = LengthCompData::id_g++;
        LengthCompData::initialized_models[this->id] = this;

    }

    virtual ~LengthCompData() {
    }

    MASObjectType sub_model_type = DATAOBJECT;
    int id;
    Rcpp::NumericVector data;
    Rcpp::NumericVector sample_size;
    std::string sex;
    double missing_values = -999;

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &selex, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
    }

    static std::map<int, LengthCompData*> initialized_models;
    typedef typename std::map<int, LengthCompData*>::iterator model_iterator;
};

std::map<int, LengthCompData*> LengthCompData::initialized_models;
int LengthCompData::id_g = 1;

class Fleet : public MASSubModel {
protected:
public:

    static int id_g;

    //    int index_data;
    //    int age_comp_data;
    //    int length_comp_data;
    int index_nll_id = -999;
    int age_comp_nll_id = -999;
    int length_comp_nll_id = -999;
    int selectivity_model_id;
    int area_id;

    std::vector<triple<int, int, int> > fishing_nortality; //id, season, area
    std::vector<triple<int, int, int> > selectivity; //id, season, area
    std::vector<std::pair<SexType, int> > index_data;
    std::vector<std::pair<SexType, int> > age_comp_data;
    std::vector<std::pair<SexType, int> > length_comp_data;
public:
    int used = false;
    double catch_fraction_of_year = 1.0;
    int id;
    std::string name;

    Fleet() {

        this->id = Fleet::id_g++;
        Fleet::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~Fleet() {
    }

    void AddIndexData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->index_data.push_back(std::make_pair(st, id));
    }

    void AddAgeCompData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->age_comp_data.push_back(std::make_pair(st, id));

    }

    void AddLengthCompData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->length_comp_data.push_back(std::make_pair(st, id));
    }

    void SetIndexCompNllComponent(int model_id) {

        this->index_nll_id = model_id;
    }

    void SetAgeCompNllComponent(int model_id) {

        this->age_comp_nll_id = model_id;
    }

    void SetLengthCompNllComponent(int model_id) {

        this->length_comp_nll_id = model_id;
    }

    void AddFishingMortality(int id, int season, int area) {

        this->fishing_nortality.push_back(
                triple<int, int, int>(id, season, area));
    }

    void AddSelectivity(int id, int season, int area) {

        this->selectivity.push_back(triple<int, int, int>(id, season, area));
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        if (this->used) {

            atl::intrusive_ptr<mas::Fleet<double> > fleet = new mas::Fleet<
                    double>();
            mas::Fleet<double> *f = fleet.get();
            f->id = this->id;
            f->catch_fraction_of_year = this->catch_fraction_of_year;
            double fraction = 0.5; //what's this?
            for (int i = 0; i < this->selectivity.size(); i++) {
                int sseason = this->selectivity[i].second;
                int sid = this->selectivity[i].first;
                int sarea = this->selectivity[i].third;
                f->operational_areas.insert(sarea);
                f->season_area_selectivity_ids[sseason][sarea] = sid;
                f->area_season_selectivity_ids[sarea][sseason] = sid;
                f->area_season_catch_fraction[sarea][sseason] = fraction;
                f->season_area_catch_fraction[sseason][sarea] = fraction;
            }

            for (int i = 0; i < this->fishing_nortality.size(); i++) {
                int sseason = this->fishing_nortality[i].second;
                int sid = this->fishing_nortality[i].first;
                int sarea = this->fishing_nortality[i].third;
                f->season_area_fishing_mortality_ids[sseason][sarea] = sid;
                f->area_season_fishing_mortality_ids[sarea][sseason] = sid;
            }

            if (this->index_nll_id != -999) {

                bool nll_exists = false;
                for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
                    if (NLLBase::nll_submodels[i]->id_ == this->index_nll_id) {
                        nll_exists = true;
                    }
                }
                if (!nll_exists) {
                    std::cout << "MAS Error: Index NLL not found!";
                    info.valid_configuration = false;
                }

                for (int i = 0; i < this->index_data.size(); i++) {
                    IndexData::model_iterator it;
                    it = IndexData::initialized_models.find(
                            this->index_data[i].second);
                    if (it != IndexData::initialized_models.end()) {
                        IndexData *data = (*it).second;
                        atl::intrusive_ptr<mas::DataObject<double> > dd =
                                new mas::DataObject<double>();
                        mas::DataObject<double> *d = dd.get();
                        d->missing_value = data->missing_values;
                        d->dimensions = 2;
                        switch (this->index_data[i].first) {
                            case MALE:
                                d->sex_type = mas::MALE;
                                break;
                            case FEMALE:
                                d->sex_type = mas::FEMALE;
                                break;
                            case UNDIFFERENTIATED:
                                d->sex_type = mas::UNDIFFERENTIATED;
                                break;
                        }

                        if (data->data.size() != info.nyears * info.nseasons
                                || data->error.size()
                                != info.nyears * info.nseasons) {
                            std::cout
                                    << "MAS Error: Index data or error vector not equal to (nseasons*nyears)\n";
                            info.valid_configuration = false;
                        } else {
                            d->imax = info.nyears;
                            d->jmax = info.nseasons;

                            for (int i = 0; i < data->data.size(); i++) {
                                d->data.push_back(data->data[i]);
                                d->observation_error.push_back(data->error[i]);
                            }
                            d->id = this->id;
                            d->name = "Index Data";
                            if (data->is_abundance) {
                                f->fishery_abundance_likelihood_component_id =
                                        this->index_nll_id;
                                d->type = mas::CATCH_ABUNDANCE;
                            } else {
                                d->type = mas::CATCH_BIOMASS;
                                f->fishery_biomass_likelihood_component_id =
                                        this->index_nll_id;
                            }
                            d->Validate();
                            info.data_dictionary[d->id].push_back(dd);
                            info.data.push_back(dd);
                        }
                    }
                }

            }

            if (this->age_comp_nll_id != -999) {
                f->fishery_age_comp_likelihood_component_id =
                        this->age_comp_nll_id;

                bool nll_exists = false;
                for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
                    if (NLLBase::nll_submodels[i]->id_
                            == this->age_comp_nll_id) {
                        nll_exists = true;
                    }
                }
                if (!nll_exists) {
                    std::cout << "MAS Error: Age Comp NLL not found!";
                    info.valid_configuration = false;
                }
                for (int i = 0; i < this->age_comp_data.size(); i++) {
                    AgeCompData::model_iterator it;
                    it = AgeCompData::initialized_models.find(
                            this->age_comp_data[i].second);
                    if (it != AgeCompData::initialized_models.end()) {
                        AgeCompData *data = (*it).second;
                        atl::intrusive_ptr<mas::DataObject<double> > dd =
                                new mas::DataObject<double>();
                        mas::DataObject<double> *d = dd.get();
                        d->missing_value = data->missing_values;
                        d->dimensions = 3;
                        switch (this->age_comp_data[i].first) {
                            case MALE:
                                d->sex_type = mas::MALE;
                                break;
                            case FEMALE:
                                d->sex_type = mas::FEMALE;
                                break;
                            case UNDIFFERENTIATED:
                                d->sex_type = mas::UNDIFFERENTIATED;
                                break;
                        }

                        if (data->data.size()
                                != info.nyears * info.nseasons
                                * info.ages.size()) {
                            std::cout
                                    << "MAS Error: Fleet Age Comp data vector not equal to (nseasons*nyears*nages) \n";
                            info.valid_configuration = false;
                        } else if (data->sample_size.size()
                                != info.nyears * info.nseasons) {
                            std::cout
                                    << "MAS Error: Fleet Age Comp error vector not equal to (nseasons*nyears) \n";
                            info.valid_configuration = false;
                        } else {

                            d->imax = info.nyears;
                            d->jmax = info.nseasons;
                            d->kmax = info.ages.size();
                            for (int i = 0; i < data->data.size(); i++) {
                                d->data.push_back(data->data[i]);
                            }
                            for (int i = 0; i < data->sample_size.size(); i++) {
                                d->sample_size.push_back(data->sample_size[i]);
                            }
                            d->id = this->id;
                            d->name = "AgeComp Data";
                            d->type = mas::CATCH_PROPORTION_AT_AGE;

                            d->Validate();
                            info.data_dictionary[d->id].push_back(dd);
                            info.data.push_back(dd);
                        }
                    }
                }

            }
            if (this->length_comp_nll_id != -999) {
                f->fishery_length_comp_likelihood_component_id =
                        this->length_comp_nll_id;

                bool nll_exists = false;
                for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
                    if (NLLBase::nll_submodels[i]->id_
                            == this->length_comp_nll_id) {
                        nll_exists = true;
                    }
                }
                if (!nll_exists) {
                    std::cout << "MAS Error: Length Comp NLL not found!";
                    info.valid_configuration = false;
                }
                for (int i = 0; i < this->length_comp_data.size(); i++) {
                    LengthCompData::model_iterator it;
                    it = LengthCompData::initialized_models.find(
                            this->length_comp_data[i].second);
                    if (it != LengthCompData::initialized_models.end()) {
                        LengthCompData *data = (*it).second;
                        atl::intrusive_ptr<mas::DataObject<double> > dd =
                                new mas::DataObject<double>();
                        mas::DataObject<double> *d = dd.get();
                        d->missing_value = data->missing_values;
                        d->dimensions = 3;
                        switch (this->length_comp_data[i].first) {
                            case MALE:
                                d->sex_type = mas::MALE;
                                break;
                            case FEMALE:
                                d->sex_type = mas::FEMALE;
                                break;
                            case UNDIFFERENTIATED:
                                d->sex_type = mas::UNDIFFERENTIATED;
                                break;
                        }

                        if (data->data.size()
                                != info.nyears * info.nseasons
                                * info.ages.size()) {
                            std::cout
                                    << "MAS Error: Fleet Length Comp data vector not equal to (nseasons*nyears*nages) \n";
                            info.valid_configuration = false;
                        } else if (data->sample_size.size()
                                != info.nyears * info.nseasons) {
                            std::cout
                                    << "MAS Error: Fleet Length Comp error vector not equal to (nseasons*nyears) \n";
                            info.valid_configuration = false;
                        } else {

                            d->imax = info.nyears;
                            d->jmax = info.nseasons;
                            d->kmax = info.ages.size();
                            for (int i = 0; i < data->data.size(); i++) {
                                d->data.push_back(data->data[i]);
                            }
                            for (int i = 0; i < data->sample_size.size(); i++) {
                                d->sample_size.push_back(data->sample_size[i]);
                            }
                            d->id = this->id;
                            d->name = "LengthComp Data";
                            d->type = mas::CATCH_MEAN_SIZE_AT_AGE;

                            d->Validate();
                            info.data_dictionary[d->id].push_back(dd);
                            info.data.push_back(dd);
                        }
                    }
                }

            }

            info.fleets[f->id] = fleet;
        } else {

            std::cout << "MAS Warning: Fleet_" << this->id
                    << " defined, but not used in the Model.\n";
            mas::mas_log << "MAS Warning: Fleet_" << this->id
                    << " defined, but not used in the Model.\n";
        }

    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void DataToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {

        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();

        for (int i = 0; i < this->index_data.size(); i++) {
            rapidjson::Value index_data(rapidjson::kObjectType);
            index_data.AddMember("data_object_type", "catch_biomass",
                    allocator);
            index_data.AddMember("name", "catch_biomass", allocator);
            index_data.AddMember("id", this->id, allocator);
            index_data.AddMember("units", "MT", allocator);
            switch (this->index_data[i].first) {
                case FEMALE:
                    index_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    index_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    index_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value error_vals(rapidjson::kArrayType);
            IndexData::model_iterator dit = IndexData::initialized_models.find(
                    this->index_data[i].second);

            if (dit != IndexData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                index_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &error = (*dit).second->error;
                MASSubModel::GenerateArrayObject(document, error_vals, error, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("observation_error", error_vals,
                        allocator);
                document.AddMember("data_object", index_data, allocator);
            } else {
                std::cout << "MAS Warning: Unable to locate index data \""
                        << this->index_data[i].second << "\".\n";
            }

        }

        rapidjson::Value age_comp_data(rapidjson::kObjectType);
        for (int i = 0; i < this->age_comp_data.size(); i++) {
            rapidjson::Value age_comp_data(rapidjson::kObjectType);
            age_comp_data.AddMember("data_object_type",
                    "catch_proportion_at_age", allocator);
            age_comp_data.AddMember("name", "catch_proportion_at_age",
                    allocator);
            age_comp_data.AddMember("id", this->id, allocator);
            age_comp_data.AddMember("units", "NA", allocator);
            switch (this->age_comp_data[i].first) {
                case FEMALE:
                    age_comp_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    age_comp_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    age_comp_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value sample_size_vals(rapidjson::kArrayType);
            AgeCompData::model_iterator dit =
                    AgeCompData::initialized_models.find(
                    this->age_comp_data[i].second);
            if (dit != AgeCompData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                age_comp_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 3,
                        nyears, nseasons, nages);
                age_comp_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &sample_size = (*dit).second->sample_size;
                MASSubModel::GenerateArrayObject(document, sample_size_vals,
                        sample_size, 2, nyears, nseasons, nages);
                age_comp_data.AddMember("sample_size", sample_size_vals,
                        allocator);
                document.AddMember("data_object", age_comp_data, allocator);
            } else {
                std::cout
                        << "MAS Warning: Unable to locate age composition data \""
                        << this->age_comp_data[i].second << "\".\n";
            }

        }
    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value fleet(rapidjson::kObjectType);
        fleet.AddMember("id", this->id, allocator);
        rapidjson::Value fishing_mortality(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->fishing_nortality.size(); i++) {
            rapidjson::Value fm_entry(rapidjson::kObjectType);
            fm_entry.AddMember("id", this->fishing_nortality[i].first,
                    allocator);
            fm_entry.AddMember("season", this->fishing_nortality[i].second,
                    allocator);
            fm_entry.AddMember("area", this->fishing_nortality[i].third,
                    allocator);
            fishing_mortality.PushBack(fm_entry, allocator);
        }

        fleet.AddMember("fishing_mortality", fishing_mortality, allocator);

        rapidjson::Value selectivity(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->selectivity.size(); i++) {
            rapidjson::Value selex_entry(rapidjson::kObjectType);
            selex_entry.AddMember("id", this->selectivity[i].first, allocator);
            selex_entry.AddMember("season", this->selectivity[i].second,
                    allocator);
            selex_entry.AddMember("area", this->selectivity[i].third,
                    allocator);
            selectivity.PushBack(selex_entry, allocator);
        }

        fleet.AddMember("selectivity", selectivity, allocator);

        rapidjson::Value likelihood(rapidjson::kArrayType);
        rapidjson::Value index(rapidjson::kObjectType);
        index.AddMember("id", this->index_nll_id, allocator);
        index.AddMember("component", "biomass_comp", allocator);
        likelihood.PushBack(index, allocator);

        rapidjson::Value age_comp(rapidjson::kObjectType);
        age_comp.AddMember("id", this->age_comp_nll_id, allocator);
        age_comp.AddMember("component", "age_comp", allocator);
        likelihood.PushBack(age_comp, allocator);

        fleet.AddMember("likelihood_components", likelihood, allocator);
        document.AddMember("fleet", fleet, allocator);
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &flt, size_t nyears, size_t nseasons, size_t nages,
            size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value fleet(rapidjson::kObjectType);
        fleet.AddMember("id", this->id, allocator);
        rapidjson::Value fishing_mortality(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->fishing_nortality.size(); i++) {
            rapidjson::Value fm_entry(rapidjson::kObjectType);
            fm_entry.AddMember("id", this->fishing_nortality[i].first,
                    allocator);
            fm_entry.AddMember("season", this->fishing_nortality[i].second,
                    allocator);
            fm_entry.AddMember("area", this->fishing_nortality[i].third,
                    allocator);
            fishing_mortality.PushBack(fm_entry, allocator);
        }

        fleet.AddMember("fishing_mortality", fishing_mortality, allocator);

        rapidjson::Value selectivity(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->selectivity.size(); i++) {
            rapidjson::Value selex_entry(rapidjson::kObjectType);
            selex_entry.AddMember("id", this->selectivity[i].first, allocator);
            selex_entry.AddMember("season", this->selectivity[i].second,
                    allocator);
            selex_entry.AddMember("area", this->selectivity[i].third,
                    allocator);
            selectivity.PushBack(selex_entry, allocator);
        }

        fleet.AddMember("selectivity", selectivity, allocator);

        rapidjson::Value likelihood(rapidjson::kArrayType);
        rapidjson::Value index(rapidjson::kObjectType);
        index.AddMember("id", this->index_nll_id, allocator);
        index.AddMember("component", "biomass_comp", allocator);
        likelihood.PushBack(index, allocator);

        rapidjson::Value age_comp(rapidjson::kObjectType);
        age_comp.AddMember("id", this->age_comp_nll_id, allocator);
        age_comp.AddMember("component", "age_comp", allocator);
        likelihood.PushBack(age_comp, allocator);

        fleet.AddMember("likelihood_components", likelihood, allocator);

        rapidjson::Value jdata(rapidjson::kArrayType);

        for (int i = 0; i < this->index_data.size(); i++) {
            rapidjson::Value index_data(rapidjson::kObjectType);
            index_data.AddMember("data_object_type", "catch_biomass",
                    allocator);
            index_data.AddMember("name", "catch_biomass", allocator);
            index_data.AddMember("id", this->id, allocator);
            index_data.AddMember("units", "MT", allocator);
            switch (this->index_data[i].first) {
                case FEMALE:
                    index_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    index_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    index_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value error_vals(rapidjson::kArrayType);
            IndexData::model_iterator dit = IndexData::initialized_models.find(
                    this->index_data[i].second);

            if (dit != IndexData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                index_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &error = (*dit).second->error;
                MASSubModel::GenerateArrayObject(document, error_vals, error, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("observation_error", error_vals,
                        allocator);
                jdata.PushBack(index_data, allocator);
            } else {
                std::cout << "MAS Warning: Unable to locate index data \""
                        << this->index_data[i].second << "\".\n";
            }

        }

        rapidjson::Value age_comp_data(rapidjson::kObjectType);
        for (int i = 0; i < this->age_comp_data.size(); i++) {
            rapidjson::Value age_comp_data(rapidjson::kObjectType);
            age_comp_data.AddMember("data_object_type",
                    "catch_proportion_at_age", allocator);
            age_comp_data.AddMember("name", "catch_proportion_at_age",
                    allocator);
            age_comp_data.AddMember("id", this->id, allocator);
            age_comp_data.AddMember("units", "NA", allocator);
            switch (this->age_comp_data[i].first) {
                case FEMALE:
                    age_comp_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    age_comp_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    age_comp_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value sample_size_vals(rapidjson::kArrayType);
            AgeCompData::model_iterator dit =
                    AgeCompData::initialized_models.find(
                    this->age_comp_data[i].second);
            if (dit != AgeCompData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                age_comp_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 3,
                        nyears, nseasons, nages);
                age_comp_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &sample_size = (*dit).second->sample_size;
                MASSubModel::GenerateArrayObject(document, sample_size_vals,
                        sample_size, 2, nyears, nseasons, nages);
                age_comp_data.AddMember("sample_size", sample_size_vals,
                        allocator);
                jdata.PushBack(age_comp_data, allocator);
            } else {
                std::cout
                        << "MAS Warning: Unable to locate age composition data \""
                        << this->age_comp_data[i].second << "\".\n";
            }

        }

        fleet.AddMember("data", jdata, allocator);

        flt.PushBack(fleet, allocator);
    }

    static std::map<int, Fleet*> initialized_models;
    typedef typename std::map<int, Fleet*>::iterator model_iterator;
};

std::map<int, Fleet*> Fleet::initialized_models;
int Fleet::id_g = 1;

class Survey : public MASSubModel {
public:

    enum Sex {
        MALE = 0, FEMALE, UNDIFFERENTIATED
    };

    static int id_g;

    std::vector<std::pair<SexType, int> > index_data;
    std::vector<std::pair<SexType, int> > age_comp_data;
    std::vector<std::pair<SexType, int> > length_comp_data;
    int index_nll_id = -999;
    int age_comp_nll_id = -999;
    int length_comp_nll_id = -999;

    std::vector<triple<int, int, int> > selectivity; //id, season, area
public:
    bool used = false;
    int id;
    std::string name;
    Parameter q;
    double survey_fraction_of_year = 1.0;

    Survey() {

        this->id = Survey::id_g++;
        Survey::initialized_models[this->id] = this;
        MASSubModel::submodels.push_back(this);
    }

    virtual ~Survey() {
    }

    void AddIndexData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->index_data.push_back(std::make_pair(st, id));
    }

    void AddAgeCompData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->age_comp_data.push_back(std::make_pair(st, id));

    }

    void AddLengthCompData(int id, std::string sex) {
        std::locale loc;

        for (std::string::size_type i = 0; i < sex.length(); ++i)
            std::tolower(sex[i], loc);

        SexType st = GetSexType(sex);
        this->length_comp_data.push_back(std::make_pair(st, id));
    }

    void SetIndexCompNllComponent(int model_id) {

        this->index_nll_id = model_id;
    }

    void SetAgeCompNllComponent(int model_id) {

        this->age_comp_nll_id = model_id;
    }

    void SetLengthCompNllComponent(int model_id) {

        this->length_comp_nll_id = model_id;
    }

    void AddSelectivity(int id, int season, int area) {

        this->selectivity.push_back(triple<int, int, int>(id, season, area));
    }

    virtual void AddToMAS(mas::Information<double> &info) {
        std::cout << "Pushing survey " << this->id
                << " to MAS model engine!!!\n\n";
        if (this->used) {
            atl::intrusive_ptr<mas::Survey<double> > survey = new mas::Survey<
                    double>();
            mas::Survey<double> *s = survey.get();
            s->id = this->id;
            s->survey_fraction_of_year = this->survey_fraction_of_year;
            double fraction = 0.5; //what's this?
            for (int i = 0; i < this->selectivity.size(); i++) {
                int sseason = this->selectivity[i].second;
                int sid = this->selectivity[i].first;
                int sarea = this->selectivity[i].third;
                s->season_area_selectivity_ids[sseason][sarea] = sid;
                s->area_season_selectivity_ids[sarea][sseason] = sid;
                s->area_season_survey_fraction[sarea][sseason] = fraction;
                s->season_area_survey_fraction[sseason][sarea] = fraction;
            }

            mas::VariableTrait<double>::SetValue(s->q, q.value);
            if (q.estimated) {

                std::stringstream ss;
                ss << "q_" << this->id;
                mas::VariableTrait<double>::SetName(s->q, ss.str());
                if (q.min != std::numeric_limits<double>::min()) {
                    mas::VariableTrait<double>::SetMinBoundary(s->q, q.min);
                }

                if (q.max != std::numeric_limits<double>::max()) {
                    mas::VariableTrait<double>::SetMaxBoundary(s->q, q.max);
                }

                survey->Register(s->q, q.phase);
                std::cout << "Survey " << survey->id
                        << " \"q\" will be estimated in phase " << q.phase
                        << "!\n";

            } else {
                std::cout << "Survey " << survey->id
                        << " \"q\" not estimated!\n";
            }

            if (this->index_nll_id != -999) {
                s->survey_biomass_likelihood_component_id = this->index_nll_id;

                bool nll_exists = false;
                for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
                    if (NLLBase::nll_submodels[i]->id_ == this->index_nll_id) {
                        nll_exists = true;
                    }
                }
                if (!nll_exists) {
                    std::cout << "MAS Error: Index NLL not found!";
                    info.valid_configuration = false;
                }

                for (int i = 0; i < this->index_data.size(); i++) {
                    IndexData::model_iterator it;
                    it = IndexData::initialized_models.find(
                            this->index_data[i].second);
                    if (it != IndexData::initialized_models.end()) {
                        IndexData *data = (*it).second;
                        atl::intrusive_ptr<mas::DataObject<double> > dd =
                                new mas::DataObject<double>();
                        mas::DataObject<double> *d = dd.get();
                        d->missing_value = data->missing_values;
                        d->dimensions = 2;
                        switch (this->index_data[i].first) {
                            case MALE:
                                d->sex_type = mas::MALE;
                                break;
                            case FEMALE:
                                d->sex_type = mas::FEMALE;
                                break;
                            case UNDIFFERENTIATED:
                                d->sex_type = mas::UNDIFFERENTIATED;
                                break;
                        }

                        if (data->data.size() != info.nyears * info.nseasons
                                || data->error.size()
                                != info.nyears * info.nseasons) {
                            std::cout
                                    << "MAS Error: Index data or error vector not equal to (nseasons*nyears)\n";
                            info.valid_configuration = false;
                        } else {
                            d->imax = info.nyears;
                            d->jmax = info.nseasons;

                            for (int i = 0; i < data->data.size(); i++) {
                                d->data.push_back(data->data[i]);
                                d->observation_error.push_back(data->error[i]);
                            }
                            d->id = this->id;
                            d->name = "Index Data";
                            if (data->is_abundance) {
                                d->type = mas::SURVEY_ABUNDANCE;
                            } else {
                                d->type = mas::SURVEY_BIOMASS;
                            }
                            d->Validate();
                            info.data_dictionary[d->id].push_back(dd);
                            info.data.push_back(dd);
                        }
                    }
                }

            }

            if (this->age_comp_nll_id != -999) {
                s->survey_age_comp_likelihood_component_id =
                        this->age_comp_nll_id;

                bool nll_exists = false;
                for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
                    if (NLLBase::nll_submodels[i]->id_
                            == this->age_comp_nll_id) {
                        nll_exists = true;
                    }
                }
                if (!nll_exists) {
                    std::cout << "MAS Error: Age Comp NLL not found!";
                    info.valid_configuration = false;
                }
                for (int i = 0; i < this->age_comp_data.size(); i++) {
                    AgeCompData::model_iterator it;
                    it = AgeCompData::initialized_models.find(
                            this->age_comp_data[i].second);
                    if (it != AgeCompData::initialized_models.end()) {
                        AgeCompData *data = (*it).second;
                        atl::intrusive_ptr<mas::DataObject<double> > dd =
                                new mas::DataObject<double>();
                        mas::DataObject<double> *d = dd.get();
                        d->missing_value = data->missing_values;
                        d->dimensions = 3;
                        switch (this->age_comp_data[i].first) {
                            case MALE:
                                d->sex_type = mas::MALE;
                                break;
                            case FEMALE:
                                d->sex_type = mas::FEMALE;
                                break;
                            case UNDIFFERENTIATED:
                                d->sex_type = mas::UNDIFFERENTIATED;
                                break;
                        }

                        if (data->data.size()
                                != info.nyears * info.nseasons
                                * info.ages.size()) {
                            std::cout
                                    << "MAS Error: Survey Age Comp data vector not equal to (nseasons*nyears*nages) \n";
                            info.valid_configuration = false;
                        } else if (data->sample_size.size()
                                != info.nyears * info.nseasons) {
                            std::cout
                                    << "MAS Error: Survey Age Comp error vector not equal to (nseasons*nyears) \n";
                            info.valid_configuration = false;
                        } else {
                            d->imax = info.nyears;
                            d->jmax = info.nseasons;
                            d->kmax = info.ages.size();
                            for (int i = 0; i < data->data.size(); i++) {
                                d->data.push_back(data->data[i]);
                            }
                            for (int i = 0; i < data->sample_size.size(); i++) {
                                d->sample_size.push_back(data->sample_size[i]);
                            }
                            d->id = this->id;
                            d->name = "AgeComp Data";
                            d->type = mas::SURVEY_PROPORTION_AT_AGE;

                            d->Validate();
                            info.data_dictionary[d->id].push_back(dd);
                            info.data.push_back(dd);
                        }
                    }
                }

            }

            info.survey_models[survey->id] = survey;
        } else {

            std::cout << "MAS Warning: Survey_" << this->id
                    << " defined, but not used in the Model.\n";
            mas::mas_log << "MAS Warning: Survey_" << this->id
                    << " defined, but not used in the Model.\n";
        }
    }

    void ExtractFromMAS(mas::Information<double> &info) {

    }

    virtual void ToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value survey(rapidjson::kObjectType);
        survey.AddMember("id", this->id, allocator);

        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value q(rapidjson::kObjectType);
        q.AddMember("value", this->q.value, allocator);
        q.AddMember("min", this->q.min, allocator);
        q.AddMember("max", this->q.max, allocator);
        q.AddMember("phase", this->q.phase, allocator);
        if (this->q.estimated) {
            q.AddMember("estimated", "true", allocator);
        } else {
            q.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("q", q, allocator);
        survey.AddMember("parameters", parameters, allocator);

        rapidjson::Value selectivity(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->selectivity.size(); i++) {
            rapidjson::Value selex_entry(rapidjson::kObjectType);
            selex_entry.AddMember("id", this->selectivity[i].first, allocator);
            selex_entry.AddMember("season", this->selectivity[i].second,
                    allocator);
            selex_entry.AddMember("area", this->selectivity[i].third,
                    allocator);
            selectivity.PushBack(selex_entry, allocator);
        }

        survey.AddMember("selectivity", selectivity, allocator);

        rapidjson::Value likelihood(rapidjson::kArrayType);
        rapidjson::Value index(rapidjson::kObjectType);
        index.AddMember("id", this->index_nll_id, allocator);
        index.AddMember("component", "biomass_comp", allocator);
        likelihood.PushBack(index, allocator);

        rapidjson::Value age_comp(rapidjson::kObjectType);
        age_comp.AddMember("id", this->age_comp_nll_id, allocator);
        age_comp.AddMember("component", "age_comp", allocator);
        likelihood.PushBack(age_comp, allocator);

        survey.AddMember("likelihood_components", likelihood, allocator);
        document.AddMember("survey", survey, allocator);
    }

    virtual void DataToJSON(rapidjson::Document &document, size_t nyears,
            size_t nseasons, size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();

        for (int i = 0; i < this->index_data.size(); i++) {
            rapidjson::Value index_data(rapidjson::kObjectType);
            index_data.AddMember("data_object_type", "survey_biomass",
                    allocator);
            index_data.AddMember("name", "survey_biomass", allocator);
            index_data.AddMember("id", this->id, allocator);
            index_data.AddMember("units", "MT", allocator);
            switch (this->index_data[i].first) {
                case FEMALE:
                    index_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    index_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    index_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value error_vals(rapidjson::kArrayType);
            IndexData::model_iterator dit = IndexData::initialized_models.find(
                    this->index_data[i].second);

            if (dit != IndexData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                index_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &error = (*dit).second->error;
                MASSubModel::GenerateArrayObject(document, error_vals, error, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("observation_error", error_vals,
                        allocator);
                document.AddMember("data_object", index_data, allocator);
            } else {
                std::cout << "MAS Warning: Unable to locate index data \""
                        << this->index_data[i].second << "\".\n";
            }

        }

        rapidjson::Value age_comp_data(rapidjson::kObjectType);
        for (int i = 0; i < this->age_comp_data.size(); i++) {
            rapidjson::Value age_comp_data(rapidjson::kObjectType);
            age_comp_data.AddMember("data_object_type",
                    "survey_proportion_at_age", allocator);
            age_comp_data.AddMember("name", "survey_proportion_at_age",
                    allocator);
            age_comp_data.AddMember("id", this->id, allocator);
            age_comp_data.AddMember("units", "NA", allocator);
            switch (this->age_comp_data[i].first) {
                case FEMALE:
                    age_comp_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    age_comp_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    age_comp_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value sample_size_vals(rapidjson::kArrayType);
            AgeCompData::model_iterator dit =
                    AgeCompData::initialized_models.find(
                    this->age_comp_data[i].second);
            if (dit != AgeCompData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                age_comp_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 3,
                        nyears, nseasons, nages);
                age_comp_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &sample_size = (*dit).second->sample_size;
                MASSubModel::GenerateArrayObject(document, sample_size_vals,
                        sample_size, 2, nyears, nseasons, nages);
                age_comp_data.AddMember("sample_size", sample_size_vals,
                        allocator);
                document.AddMember("data_object", age_comp_data, allocator);
            } else {
                std::cout
                        << "MAS Warning: Unable to locate age composition data \""
                        << this->age_comp_data[i].second << "\".\n";
            }

        }
    }

    virtual void AddToEMInputs(rapidjson::Document &document,
            rapidjson::Value &srvy, size_t nyears, size_t nseasons,
            size_t nages, size_t nareas) {
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();
        rapidjson::Value survey(rapidjson::kObjectType);
        survey.AddMember("id", this->id, allocator);

        rapidjson::Value parameters(rapidjson::kObjectType);
        rapidjson::Value q(rapidjson::kObjectType);
        q.AddMember("value", this->q.value, allocator);
        q.AddMember("min", this->q.min, allocator);
        q.AddMember("max", this->q.max, allocator);
        q.AddMember("phase", this->q.phase, allocator);
        if (this->q.estimated) {
            q.AddMember("estimated", "true", allocator);
        } else {
            q.AddMember("estimated", "false", allocator);
        }
        parameters.AddMember("q", q, allocator);
        survey.AddMember("parameters", parameters, allocator);

        rapidjson::Value selectivity(rapidjson::kArrayType);
        //id, season, area
        for (int i = 0; i < this->selectivity.size(); i++) {
            rapidjson::Value selex_entry(rapidjson::kObjectType);
            selex_entry.AddMember("id", this->selectivity[i].first, allocator);
            selex_entry.AddMember("season", this->selectivity[i].second,
                    allocator);
            selex_entry.AddMember("area", this->selectivity[i].third,
                    allocator);
            selectivity.PushBack(selex_entry, allocator);
        }

        survey.AddMember("selectivity", selectivity, allocator);

        rapidjson::Value likelihood(rapidjson::kArrayType);
        rapidjson::Value index(rapidjson::kObjectType);
        index.AddMember("id", this->index_nll_id, allocator);
        index.AddMember("component", "biomass_comp", allocator);
        likelihood.PushBack(index, allocator);

        rapidjson::Value age_comp(rapidjson::kObjectType);
        age_comp.AddMember("id", this->age_comp_nll_id, allocator);
        age_comp.AddMember("component", "age_comp", allocator);
        likelihood.PushBack(age_comp, allocator);

        survey.AddMember("likelihood_components", likelihood, allocator);

        rapidjson::Value jdata(rapidjson::kArrayType);

        for (int i = 0; i < this->index_data.size(); i++) {
            rapidjson::Value index_data(rapidjson::kObjectType);
            index_data.AddMember("data_object_type", "survey_biomass",
                    allocator);
            index_data.AddMember("name", "survey_biomass", allocator);
            index_data.AddMember("id", this->id, allocator);
            index_data.AddMember("units", "MT", allocator);
            switch (this->index_data[i].first) {
                case FEMALE:
                    index_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    index_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    index_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value error_vals(rapidjson::kArrayType);
            IndexData::model_iterator dit = IndexData::initialized_models.find(
                    this->index_data[i].second);

            if (dit != IndexData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                index_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &error = (*dit).second->error;
                MASSubModel::GenerateArrayObject(document, error_vals, error, 2,
                        nyears, nseasons, nages);
                index_data.AddMember("observation_error", error_vals,
                        allocator);
                jdata.PushBack(index_data, allocator);
            } else {
                std::cout << "MAS Warning: Unable to locate index data \""
                        << this->index_data[i].second << "\".\n";
            }

        }

        rapidjson::Value age_comp_data(rapidjson::kObjectType);
        for (int i = 0; i < this->age_comp_data.size(); i++) {
            rapidjson::Value age_comp_data(rapidjson::kObjectType);
            age_comp_data.AddMember("data_object_type",
                    "survey_proportion_at_age", allocator);
            age_comp_data.AddMember("name", "survey_proportion_at_age",
                    allocator);
            age_comp_data.AddMember("id", this->id, allocator);
            age_comp_data.AddMember("units", "NA", allocator);
            switch (this->age_comp_data[i].first) {
                case FEMALE:
                    age_comp_data.AddMember("sex", "female", allocator);
                    break;
                case MALE:
                    age_comp_data.AddMember("sex", "male", allocator);
                    break;
                case UNDIFFERENTIATED:
                    age_comp_data.AddMember("sex", "undifferentiated", allocator);
                    break;
                default:
                    std::cout << "MAS Warning: Unknown sex type for index data.\n";
            }
            rapidjson::Value vals(rapidjson::kArrayType);
            rapidjson::Value sample_size_vals(rapidjson::kArrayType);
            AgeCompData::model_iterator dit =
                    AgeCompData::initialized_models.find(
                    this->age_comp_data[i].second);
            if (dit != AgeCompData::initialized_models.end()) {

                double missing_values = (*dit).second->missing_values;
                age_comp_data.AddMember("missing_values", missing_values,
                        allocator);

                Rcpp::NumericVector &data = (*dit).second->data;
                MASSubModel::GenerateArrayObject(document, vals, data, 3,
                        nyears, nseasons, nages);
                age_comp_data.AddMember("values", vals, allocator);

                Rcpp::NumericVector &sample_size = (*dit).second->sample_size;
                MASSubModel::GenerateArrayObject(document, sample_size_vals,
                        sample_size, 2, nyears, nseasons, nages);
                age_comp_data.AddMember("sample_size", sample_size_vals,
                        allocator);
                jdata.PushBack(age_comp_data, allocator);
            } else {
                std::cout
                        << "MAS Warning: Unable to locate age composition data \""
                        << this->age_comp_data[i].second << "\".\n";
            }

        }

        srvy.AddMember("data", jdata, allocator);
        srvy.PushBack(survey, allocator);
    }

    static std::map<int, Survey*> initialized_models;
    typedef typename std::map<int, Survey*>::iterator model_iterator;
};

std::map<int, Survey*> Survey::initialized_models;
int Survey::id_g = 1;

class MASModel {
    bool initialized = false;
    std::set<int> fleets;
    std::set<int> surveys;
    std::set<int> populations;


    //data produced from the operating model
    std::vector<atl::intrusive_ptr<IndexData> > om_index_data;
    std::vector<atl::intrusive_ptr<AgeCompData> > om_age_comp_data;

private:

    /**
     * Create the actual MAS model
     */
    void Initialize() {
        if (!this->initialized) {

        }
    }

    std::shared_ptr<mas::MASObjectiveFunction<double> > mas;
public:
    bool compute_variance_for_derived_quantities = true;
    int nyears;
    int nseasons;
    int nages;
    int max_line_searches = 50;
    int max_iterations = 1000;
    int print_interval = 10;
    double tolerance = 1e-4;
    double spawning_season_offset = 0.0;
    double catch_season_offset = 0.0;
    double survey_season_offset = 0.0;
    double extended_plus_group = 0;
    Rcpp::NumericVector ages;

    virtual ~MASModel() {

    }

    int GetNAges() const {

        return nages;
    }

    void SetNAges(int nages) {

        this->nages = nages;
    }

    int GetNSeasons() const {

        return nseasons;
    }

    void SetNSeasons(int nseasons) {

        this->nseasons = nseasons;
    }

    int GetNYears() const {

        return nyears;
    }

    void SetNYears(int nyears) {

        this->nyears = nyears;
    }

    void Run() {

        mas = std::make_shared<mas::MASObjectiveFunction<double> >();
        mas->compute_variance_for_derived_quantities = this->compute_variance_for_derived_quantities;

        if (this->nages == 0) {
            std::cout << "MAS error: nages = 0\n";
            return;
        }
        if (this->nyears == 0) {
            std::cout << "MAS error: nyears = 0\n";
            return;
        }

        if (this->extended_plus_group == 0) {
            std::cout << "MAS error: extended_plus_group = 0\n";
            return;
        } else {
            mas::Subpopulation<double>::length_weight_key_carryout =
                    this->extended_plus_group;
        }

        mas->mas_instance.info.nyears = this->nyears;
        mas->mas_instance.info.nseasons = this->nseasons;
        mas->mas_instance.info.spawning_season_offset =
                this->spawning_season_offset;
        mas->mas_instance.info.survey_fraction_of_year =
                this->survey_season_offset;
        mas->mas_instance.info.catch_fraction_of_year =
                this->catch_season_offset;

        typedef typename mas::VariableTrait<double>::variable variable;
        mas->mas_instance.info.ages.resize(this->nages);
        mas->mas_instance.info.ages_real.resize(this->nages);
        mas::GrowthBase<double>::ages.clear();
        mas::GrowthBase<double>::ages_to_intrpolate.clear();
        for (int i = 0; i < nages; i++) {
            mas->mas_instance.info.ages[i] = variable(this->ages[i]);
            std::cout << mas->mas_instance.info.ages[i] << "  ";
            mas->mas_instance.info.ages_real[i] = (this->ages[i]);
            mas::GrowthBase<double>::ages.push_back(this->ages[i]);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(this->ages[i]);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->catch_season_offset);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->survey_season_offset);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->spawning_season_offset);

        }

        std::cout << "\n\n";
        for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
            NLLBase::nll_submodels[i]->AddToMAS(mas->mas_instance.info);
        }

        for (int i = 0; i < MASSubModel::submodels.size(); i++) {
            MASSubModel::submodels[i]->AddToMAS(mas->mas_instance.info);
        }

        mas->mas_instance.info.CreateModel();
        mas->Initialize();

        try {

            if (mas->mas_instance.info.valid_configuration == true) {

                atl::LBFGS<double> min;
                min.SetPrintWidth(2);
                min.max_line_searches = this->max_line_searches;
                min.SetTolerance(this->tolerance);
                min.max_iterations = this->max_iterations;
                min.print_interval = this->print_interval;
                min.SetObjectiveFunction(mas.get());
                min.Run();
                mas->Finalize();
            } else {

                std::cout
                        << "MAS Error: Invalid Model Configuration, see mas.log\n";
            }
        } catch (...) {

        }

    }

    void RunOM() {

        mas = std::make_shared<mas::MASObjectiveFunction<double> >();
        if (this->nages == 0) {
            std::cout << "MAS error: nages = 0\n";
            return;
        }
        if (this->nyears == 0) {
            std::cout << "MAS error: nyears = 0\n";
            return;
        }

        if (this->extended_plus_group == 0) {
            std::cout << "MAS error: extended_plus_group = 0\n";
            return;
        } else {
            mas::Subpopulation<double>::length_weight_key_carryout =
                    this->extended_plus_group;
        }

        mas->mas_instance.info.nyears = this->nyears;
        mas->mas_instance.info.nseasons = this->nseasons;
        std::cout << "integer check: " << mas->mas_instance.info.nyears << "  "
                << mas->mas_instance.info.nseasons << "\n";
        mas->mas_instance.info.spawning_season_offset =
                this->spawning_season_offset;
        mas->mas_instance.info.survey_fraction_of_year =
                this->survey_season_offset;
        mas->mas_instance.info.catch_fraction_of_year =
                this->catch_season_offset;

        typedef typename mas::VariableTrait<double>::variable variable;
        mas->mas_instance.info.ages.resize(this->nages);
        mas->mas_instance.info.ages_real.resize(this->nages);
        mas::GrowthBase<double>::ages.clear();
        mas::GrowthBase<double>::ages_to_intrpolate.clear();
        for (int i = 0; i < nages; i++) {
            mas->mas_instance.info.ages[i] = variable(this->ages[i]);

            mas->mas_instance.info.ages_real[i] = (this->ages[i]);
            mas::GrowthBase<double>::ages.push_back(this->ages[i]);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(this->ages[i]);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->catch_season_offset);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->survey_season_offset);
            mas::GrowthBase<double>::ages_to_intrpolate.insert(
                    this->ages[i] + this->spawning_season_offset);

        }

        for (int i = 0; i < NLLBase::nll_submodels.size(); i++) {
            NLLBase::nll_submodels[i]->AddToMAS(mas->mas_instance.info);
        }

        for (int i = 0; i < MASSubModel::submodels.size(); i++) {
            MASSubModel::submodels[i]->AddToMAS(mas->mas_instance.info);
        }

        mas->mas_instance.info.CreateModel();
        mas->Initialize();

        mas->mas_instance.RunOperationalModel();

        mas->Finalize();

        //transfer derived values from MAS to RMAS
        Fleet::model_iterator it;

        for (it = Fleet::initialized_models.begin();
                it != Fleet::initialized_models.end(); ++it) {
            Fleet *f = (*it).second;

            atl::intrusive_ptr<IndexData> fleet_index_data;
            atl::intrusive_ptr<AgeCompData> fleet_age_comp_data;

            int id = f->id;
            atl::intrusive_ptr<mas::DataObject<double> > data =
                    mas->mas_instance.info.fleets[id]->catch_biomass_data;

            fleet_index_data = new IndexData();
            //            fleet_index_data->id = data->id;
            fleet_index_data->data = data->data;
            fleet_index_data->error = data->observation_error;
            fleet_index_data->sex = "undifferentiated";
            this->om_index_data.push_back(fleet_index_data);
            //            IndexData::initialized_models[fleet_index_data->id] = fleet_index_data.get();
            //add back to initialized list
            f->AddIndexData(fleet_index_data->id, "undifferentiated");

            atl::intrusive_ptr<mas::DataObject<double> > data2 =
                    mas->mas_instance.info.fleets[id]->catch_proportion_at_age_data;

            fleet_age_comp_data = new AgeCompData();
            //            fleet_age_comp_data->id = data2->id;
            fleet_age_comp_data->data = data2->data;
            fleet_age_comp_data->sample_size = data2->sample_size;
            fleet_age_comp_data->sex = "undifferentiated";
            this->om_age_comp_data.push_back(fleet_age_comp_data);
            //            AgeCompData::initialized_models[fleet_age_comp_data->id] = fleet_age_comp_data.get();
            //add back to initialized list
            f->AddAgeCompData(fleet_age_comp_data->id, "undifferentiated");
        }

        Survey::model_iterator sit;
        for (sit = Survey::initialized_models.begin();
                sit != Survey::initialized_models.end(); ++sit) {

            std::cout << "here" << std::endl;

            Survey *s = (*sit).second;

            atl::intrusive_ptr<IndexData> survey_index_data;
            atl::intrusive_ptr<AgeCompData> survey_age_comp_data;

            int id = s->id;
            atl::intrusive_ptr<mas::DataObject<double> > data =
                    mas->mas_instance.info.survey_models[id]->survey_biomass_data;

            survey_index_data = new IndexData();
            //            survey_index_data->id = data->id;
            survey_index_data->data = data->data;
            survey_index_data->error = data->observation_error;
            survey_index_data->sex = "undifferentiated";
            this->om_index_data.push_back(survey_index_data);
            //            IndexData::initialized_models[survey_index_data->id] = survey_index_data.get();
            //add back to initialized list
            s->AddIndexData(survey_index_data->id, "undifferentiated");

            atl::intrusive_ptr<mas::DataObject<double> > data2 =
                    mas->mas_instance.info.survey_models[id]->survey_proportion_at_age_data;

            survey_age_comp_data = new AgeCompData();
            //            survey_age_comp_data->id = data2->id;
            survey_age_comp_data->data = data2->data;
            survey_age_comp_data->sample_size = data2->sample_size;
            survey_age_comp_data->sex = "undifferentiated";
            this->om_age_comp_data.push_back(survey_age_comp_data);
            //            AgeCompData::initialized_models[survey_age_comp_data->id] = survey_age_comp_data.get();
            //add back to initialized list
            s->AddAgeCompData(survey_age_comp_data->id, "undifferentiated");
        }

    }

    void AddFleet(int id) {
        Fleet::model_iterator it = Fleet::initialized_models.find(id);
        if (it != Fleet::initialized_models.end()) {
            (

                    *it).second->used = true;
            this->fleets.insert(id);
        }

    }

    void AddSurvey(int id) {
        Survey::model_iterator it = Survey::initialized_models.find(id);
        if (it != Survey::initialized_models.end()) {
            (

                    *it).second->used = true;
            this->surveys.insert(id);
        }
    }

    void AddPopulation(int id) {

        this->populations.insert(id);
    }

    std::string GetOuptput() {

        mas::JSONOutputGenerator<double> json;
        mas->SetVarianceCovariance();

        return json.GenerateOutput(mas->mas_instance);
    }

    std::string GetJSONData() {
        rapidjson::Document document;
        document.SetObject();
        int nareas = Area::initialized_models.size();
        Fleet::model_iterator fit;
        for (fit = Fleet::initialized_models.begin();
                fit != Fleet::initialized_models.end(); ++fit) {
            (*fit).second->DataToJSON(document, nyears, nseasons, nages,
                    nareas);
        }

        Survey::model_iterator sit;
        for (sit = Survey::initialized_models.begin();
                sit != Survey::initialized_models.end(); ++sit) {
            (*sit).second->DataToJSON(document, nyears, nseasons, nages,
                    nareas);
        }

        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        document.Accept(writer);
        return buffer.GetString();

    }

    std::string GetJSONConfig() {
        rapidjson::Document document;
        document.SetObject();
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();

        document.AddMember("years", this->nyears, allocator);
        document.AddMember("seasons", this->nseasons, allocator);
        document.AddMember("spawning_season_offset",
                this->spawning_season_offset, allocator);
        document.AddMember("catch_season_offset", this->catch_season_offset,
                allocator);
        document.AddMember("survey_season_offset", this->survey_season_offset,
                allocator);
        document.AddMember("extended_plus_group",
                static_cast<int> (this->extended_plus_group), allocator);

        rapidjson::Value agev(rapidjson::kArrayType);
        for (int i = 0; i < this->ages.size(); i++) {
            agev.PushBack(this->ages[i], allocator);
        }
        document.AddMember("ages", agev, allocator);

        size_t nareas = Area::initialized_models.size();
        for (int i = 0; i < MASSubModel::submodels.size(); i++) {
            if (MASSubModel::submodels[i]->sub_model_type != DATAOBJECT) {
                MASSubModel::submodels[i]->ToJSON(document, nyears, nseasons,
                        nages, nareas);
            }
        }

        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        document.Accept(writer);
        std::cout << buffer.GetString();
        return buffer.GetString();

    }

    std::string GetEMInputs() {
        rapidjson::Document document;
        document.SetObject();
        rapidjson::Document::AllocatorType &allocator = document.GetAllocator();

        document.AddMember("years", this->nyears, allocator);
        document.AddMember("seasons", this->nseasons, allocator);
        document.AddMember("spawning_season_offset",
                this->spawning_season_offset, allocator);
        document.AddMember("catch_season_offset", this->catch_season_offset,
                allocator);
        document.AddMember("survey_season_offset", this->survey_season_offset,
                allocator);
        document.AddMember("extended_plus_group",
                static_cast<int> (this->extended_plus_group), allocator);

        rapidjson::Value agev(rapidjson::kArrayType);
        for (int i = 0; i < this->ages.size(); i++) {
            agev.PushBack(this->ages[i], allocator);
        }
        document.AddMember("ages", agev, allocator);

        size_t nareas = Area::initialized_models.size();

        rapidjson::Value areas(rapidjson::kArrayType);
        typename Area::model_iterator ait;

        for (ait = Area::initialized_models.begin();
                ait != Area::initialized_models.end(); ++ait) {
            (*ait).second->AddToEMInputs(document, areas, this->nyears,
                    this->nseasons, this->nages, nareas);
        }
        document.AddMember("areas", areas, allocator);

        //        rapidjson::Value recruitment(rapidjson::kArrayType);
        //        typename RecruitmentBase::model_iterator rit;
        //
        //        for (rit = RecruitmentBase::initialized_models.begin(); rit != RecruitmentBase::initialized_models.end(); ++rit) {
        //            (*rit).second->AddToEMInputs(document, recruitment, this->nyears, this->nseasons, this->nages, nareas);
        //        }
        //        document.AddMember("recruitment", recruitment, allocator);
        //
        //        rapidjson::Value growth(rapidjson::kArrayType);
        //
        //        typename GrowthBase::

        rapidjson::Value maturity(rapidjson::kArrayType);
        rapidjson::Value mortality(rapidjson::kArrayType);
        rapidjson::Value initial_deviations(rapidjson::kArrayType);
        rapidjson::Value movement(rapidjson::kArrayType);

        rapidjson::Value selectivity(rapidjson::kArrayType);
        rapidjson::Value fishing_mortality(rapidjson::kArrayType);
        rapidjson::Value fleets(rapidjson::kArrayType);
        rapidjson::Value surveys(rapidjson::kArrayType);

        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        document.Accept(writer);
        std::cout << buffer.GetString();
        return buffer.GetString();

    }

    Rcpp::List GetParameterEstimates() {

        Rcpp::List l;

    }

    void Reset() {
        LogisticSelectivity::id_g = 1;
        LogisticSelectivity::initialized_models.clear();

        DoubleLogisticSelectivity::id_g = 1;
        DoubleLogisticSelectivity::initialized_models.clear();

        AgeBasedSelectivity::id_g = 1;
        AgeBasedSelectivity::initialized_models.clear();

        FishingMortality::id_g = 1;
        FishingMortality::initialized_models.clear();

        NaturalMortality::id_g = 1;
        NaturalMortality::initialized_models.clear();

        InitialDeviations::id_g = 1;
        InitialDeviations::initialized_models.clear();

        RickerRecruitment::id_g = 1;
        RickerRecruitment::initialized_models.clear();

        BevertonHoltRecruitment::id_g = 1;
        BevertonHoltRecruitment::initialized_models.clear();

        VonBertalanffy::id_g = 1;
        VonBertalanffy::initialized_models.clear();

        VonBertalanffyModified::id_g = 1;
        VonBertalanffyModified::initialized_models.clear();

        Area::id_g = 1;
        Area::initialized_models.clear();

        Movement::id_g = 1;
        Movement::initialized_models.clear();

        Maturity::id_g = 1;
        Maturity::initialized_models.clear();

        Population::id_g = 1;
        Population::initialized_models.clear();

        IndexData::id_g = 1;
        IndexData::initialized_models.clear();

        AgeCompData::id_g = 1;
        AgeCompData::initialized_models.clear();

        LengthCompData::id_g = 1;
        LengthCompData::initialized_models.clear();

        Fleet::id_g = 1;
        Fleet::initialized_models.clear();

        Survey::id_g = 1;
        Survey::initialized_models.clear();

        NLLBase::id_g = 1;
        NLLBase::nll_submodels.clear();
        MASSubModel::submodels.clear();
    }

};

RCPP_EXPOSED_CLASS(Parameter)
RCPP_EXPOSED_CLASS(LogisticSelectivity)
RCPP_EXPOSED_CLASS(DoubleLogisticSelectivity)
RCPP_EXPOSED_CLASS(FishingMortality)
RCPP_EXPOSED_CLASS(NaturalMortality)
RCPP_EXPOSED_CLASS(RickerRecruitment)
RCPP_EXPOSED_CLASS(BevertonHoltRecruitment)
RCPP_EXPOSED_CLASS(BevertonHoltRecruitmentAlt)
RCPP_EXPOSED_CLASS(VonBertalanffy)
RCPP_EXPOSED_CLASS(VonBertalanffyModified)
RCPP_EXPOSED_CLASS(Area)
RCPP_EXPOSED_CLASS(Movement)
RCPP_EXPOSED_CLASS(Maturity)
RCPP_EXPOSED_CLASS(Population)
RCPP_EXPOSED_CLASS(Lognormal)
RCPP_EXPOSED_CLASS(DirichletMultinomial)
RCPP_EXPOSED_CLASS(DirichletMultinomialRobust)
RCPP_EXPOSED_CLASS(Multinomial)
RCPP_EXPOSED_CLASS(MultinomialRobust)
RCPP_EXPOSED_CLASS(Fleet)
RCPP_EXPOSED_CLASS(Survey)
RCPP_EXPOSED_CLASS(IndexData)
RCPP_EXPOSED_CLASS(AgeCompData)
RCPP_EXPOSED_CLASS(LengthCompData)
RCPP_EXPOSED_CLASS(MASModel)

RCPP_MODULE(ht4sa) {
    class_<Parameter>("Parameter")
            .constructor()
            .constructor<double>()
            .constructor<Parameter>()
            .field("value", &Parameter::value)
            .field("min", &Parameter::min)
            .field("max", &Parameter::max)
            .field("estimated", &Parameter::estimated)
            .field("phase", &Parameter::phase)
            .field("lambda", &Parameter::lambda);

    class_<LogisticSelectivity>("LogisticSelectivity")
            .constructor()
            .field("a50", &LogisticSelectivity::a50)
            .field("slope", &LogisticSelectivity::slope)
            .field("sigma", &LogisticSelectivity::sigma)
            .field("sigma2", &LogisticSelectivity::sigma2)
            .field("cv", &LogisticSelectivity::cv)
            .field("id", &LogisticSelectivity::id)
            ;

    class_<AgeBasedSelectivity>("AgeBasedSelectivity")
            .constructor()
            .field("values", &AgeBasedSelectivity::values)
            .field("estimate_age", &AgeBasedSelectivity::estimate_age)
            .field("sigma", &AgeBasedSelectivity::sigma)
            .field("sigma2", &AgeBasedSelectivity::sigma2)
            .field("cv", &AgeBasedSelectivity::cv)
            .field("id", &AgeBasedSelectivity::id)
            .field("estimated", &AgeBasedSelectivity::estimated)
            .field("phase", &AgeBasedSelectivity::phase)
            ;

    class_<DoubleLogisticSelectivity>("DoubleLogisticSelectivity")
            .constructor()
            .field("alpha_asc", &DoubleLogisticSelectivity::alpha_asc)
            .field("beta_asc", &DoubleLogisticSelectivity::beta_asc)
            .field("alpha_desc", &DoubleLogisticSelectivity::alpha_desc)
            .field("beta_desc", &DoubleLogisticSelectivity::beta_desc)
            .field("sigma", &DoubleLogisticSelectivity::sigma)
            .field("sigma2", &DoubleLogisticSelectivity::sigma2)
            .field("cv", &DoubleLogisticSelectivity::cv)
            .field("id", &DoubleLogisticSelectivity::id);
    class_<FishingMortality>("FishingMortality")
            .constructor()
            .method("SetValues", &FishingMortality::SetValues)
            .field("id", &FishingMortality::id)
            .field("values", &FishingMortality::values)
            .field("min", &FishingMortality::min)
            .field("max", &FishingMortality::max)
            .field("estimate", &FishingMortality::estimate)
            .field("phase", &FishingMortality::phase)
            ;

    class_<NaturalMortality>("NaturalMortality")
            .constructor()
            .method("SetValues", &NaturalMortality::SetValues)
            .field("id", &NaturalMortality::id)
            .field("values", &NaturalMortality::values)
            .field("estimate", &NaturalMortality::estimate)
            .field("phase", &NaturalMortality::phase)
            ;

    class_<InitialDeviations>("InitialDeviations")
            .constructor()
            .method("SetValues", &InitialDeviations::SetValues)
            .field("id", &InitialDeviations::id)
            .field("values", &InitialDeviations::values)
            .field("estimate", &InitialDeviations::estimate)
            .field("phase", &InitialDeviations::phase)
            ;

    class_<RickerRecruitment>("RickerRecruitment")
            .constructor()
            .method("SetDeviations", &RickerRecruitment::SetDeviations)
            .field("R0", &RickerRecruitment::R0)
            .field("use_bias_correction", &RickerRecruitment::use_bias_correction)
            .method("Evaluate", &RickerRecruitment::Evaluate)
            .field("alpha", &RickerRecruitment::alpha)
            .field("constrained_deviations", &RickerRecruitment::constrained_deviations)
            .field("estimate_deviations", &RickerRecruitment::estimate_deviations)
            .field("deviation_phase", &RickerRecruitment::deviation_phase)
            .field("deviations_min", &RickerRecruitment::deviations_min)
            .field("deviations_max", &RickerRecruitment::deviations_max)
            .field("id", &RickerRecruitment::id)
            ;

    class_<BevertonHoltRecruitment>("BevertonHoltRecruitment")
            .constructor()
            .method("SetDeviations", &BevertonHoltRecruitment::SetDeviations)
            .method("Evaluate", &BevertonHoltRecruitment::Evaluate)
            .field("R0", &BevertonHoltRecruitment::R0)
            .field("use_bias_correction", &BevertonHoltRecruitment::use_bias_correction)
            .field("sigma_r", &BevertonHoltRecruitment::sigma_r)
            .field("h", &BevertonHoltRecruitment::h)
            .field("constrained_deviations", &BevertonHoltRecruitment::constrained_deviations)
            .field("estimate_deviations", &BevertonHoltRecruitment::estimate_deviations)
            .field("deviation_phase", &BevertonHoltRecruitment::deviation_phase)
            .field("deviations_min", &BevertonHoltRecruitment::deviations_min)
            .field("deviations_max", &BevertonHoltRecruitment::deviations_max)
            .field("id", &BevertonHoltRecruitment::id)
            ;

    class_<BevertonHoltRecruitmentAlt>("BevertonHoltRecruitmentAlt")
            .constructor()
            .method("SetDeviations", &BevertonHoltRecruitmentAlt::SetDeviations)
            .method("Evaluate", &BevertonHoltRecruitmentAlt::Evaluate)
            .field("R0", &BevertonHoltRecruitmentAlt::R0)
            .field("use_bias_correction", &BevertonHoltRecruitmentAlt::use_bias_correction)
            .field("sigma_r", &BevertonHoltRecruitmentAlt::sigma_r)
            .field("h", &BevertonHoltRecruitmentAlt::h)
            .field("constrained_deviations", &BevertonHoltRecruitmentAlt::constrained_deviations)
            .field("estimate_deviations", &BevertonHoltRecruitmentAlt::estimate_deviations)
            .field("deviation_phase", &BevertonHoltRecruitmentAlt::deviation_phase)
            .field("deviations_min", &BevertonHoltRecruitmentAlt::deviations_min)
            .field("deviations_max", &BevertonHoltRecruitmentAlt::deviations_max)
            .field("id", &BevertonHoltRecruitmentAlt::id)
            ;

    class_<VonBertalanffy>("VonBertalanffy")
            .constructor()
            .method("SetUndifferentiatedWeightAtSeasonStart", &VonBertalanffy::SetUndifferentiatedWeightAtSeasonStart)
            .method("SetUndifferentiatedWeightAtSpawning", &VonBertalanffy::SetUndifferentiatedWeightAtSpawning)
            .method("SetUndifferentiatedCatchWeight", &VonBertalanffy::SetUndifferentiatedCatchWeight)
            .method("SetUndifferentiatedSurveyWeight", &VonBertalanffy::SetUndifferentiatedSurveyWeight)
            .method("SetMaleWeightAtSeasonStart", &VonBertalanffy::SetMaleWeightAtSeasonStart)
            .method("SetMaleWeightAtSpawning", &VonBertalanffy::SetMaleWeightAtSpawning)
            .method("SetMaleCatchWeight", &VonBertalanffy::SetMaleCatchWeight)
            .method("SetMaleSurveyWeight", &VonBertalanffy::SetMaleSurveyWeight)
            .method("SetFemaleWeightAtSeasonStart", &VonBertalanffy::SetFemaleWeightAtSeasonStart)
            .method("SetFemaleWeightAtSpawning", &VonBertalanffy::SetFemaleWeightAtSpawning)
            .method("SetFemaleCatchWeight", &VonBertalanffy::SetFemaleCatchWeight)
            .method("SetFemaleSurveyWeight", &VonBertalanffy::SetFemaleSurveyWeight)
            .field("a_max", &VonBertalanffy::a_max)
            .field("a_min", &VonBertalanffy::a_min)
            .field("alpha_f", &VonBertalanffy::alpha_f)
            .field("alpha_m", &VonBertalanffy::alpha_m)
            .field("beta_f", &VonBertalanffy::beta_f)
            .field("beta_m", &VonBertalanffy::beta_m)
            .field("k", &VonBertalanffy::k)
            .field("l_inf", &VonBertalanffy::l_inf)
            .field("id", &VonBertalanffy::id)
            ;

    class_<VonBertalanffyModified>("VonBertalanffyModified")
            .constructor()
            .method("SetUndifferentiatedWeightAtSeasonStart", &VonBertalanffyModified::SetUndifferentiatedWeightAtSeasonStart)
            .method("SetUndifferentiatedWeightAtSpawning", &VonBertalanffyModified::SetUndifferentiatedWeightAtSpawning)
            .method("SetUndifferentiatedCatchWeight", &VonBertalanffyModified::SetUndifferentiatedCatchWeight)
            .method("SetUndifferentiatedSurveyWeight", &VonBertalanffyModified::SetUndifferentiatedSurveyWeight)
            .method("SetMaleWeightAtSeasonStart", &VonBertalanffyModified::SetMaleWeightAtSeasonStart)
            .method("SetMaleWeightAtSpawning", &VonBertalanffyModified::SetMaleWeightAtSpawning)
            .method("SetMaleCatchWeight", &VonBertalanffyModified::SetMaleCatchWeight)
            .method("SetMaleSurveyWeight", &VonBertalanffyModified::SetMaleSurveyWeight)
            .method("SetFemaleWeightAtSeasonStart", &VonBertalanffyModified::SetFemaleWeightAtSeasonStart)
            .method("SetFemaleWeightAtSpawning", &VonBertalanffyModified::SetFemaleWeightAtSpawning)
            .method("SetFemaleCatchWeight", &VonBertalanffyModified::SetFemaleCatchWeight)
            .method("SetFemaleSurveyWeight", &VonBertalanffyModified::SetFemaleSurveyWeight)
            .field("a_max", &VonBertalanffyModified::a_max)
            .field("a_min", &VonBertalanffyModified::a_min)
            .field("alpha_f", &VonBertalanffyModified::alpha_f)
            .field("alpha_m", &VonBertalanffyModified::alpha_m)
            .field("beta_f", &VonBertalanffyModified::beta_f)
            .field("beta_m", &VonBertalanffyModified::beta_m)
            .field("c", &VonBertalanffyModified::c)
            .field("lmin", &VonBertalanffyModified::lmin)
            .field("lmax", &VonBertalanffyModified::lmax)
            .field("l_inf", &VonBertalanffyModified::l_inf)
            .field("id", &VonBertalanffyModified::id)
            ;

    class_<Area>("Area")
            .constructor()
            .field("id", &Area::id)
            .field("name", &Area::name)
            ;

    class_<Movement>("Movement")
            .constructor()
            .field("id", &Movement::id)
            .field("connectivity_males", &Movement::connectivity_males)
            .field("connectivity_females", &Movement::connectivity_females)
            .field("connectivity_recruits", &Movement::connectivity_recruits)
            .field("estimate_movement", &Movement::estimate_movement)
            .field("movement_phase", &Movement::movement_phase)
            ;

    class_<Maturity>("Maturity")
            .constructor()
            .field("id", &Maturity::id)
            .field("values", &Maturity::values)
            ;

    class_<Population>("Population")
            .constructor()
            .field("id", &Population::id)
            .field("sex_ratio", &Population::sex_ratio)
            .field("spawning_season_offset", &Population::spawning_season_offset)
            .method("SetGrowth", &Population::SetGrowth)
            .method("AddMaturity", &Population::AddMaturity)
            .method("AddMovement", &Population::AddMovement)
            .method("AddNaturalMortality", &Population::AddNaturalMortality)
            .method("AddRecruitment", &Population::AddRecruitment)
            .method("SetInitialDeviations", &Population::SetInitialDeviations)
            ;
    class_<Lognormal>("Lognormal")
            .constructor()
            .field("id", &Lognormal::id)
            .field("use_bias_correction", &Lognormal::use_bias_correction)
            .method("SetLambdaValues", &Lognormal::SetLambdaValues)
            ;
    class_<DirichletMultinomial>("DirichletMultinomial")
            .constructor()
            .field("id", &DirichletMultinomial::id)
            .method("SetLambdaValues", &DirichletMultinomial::SetLambdaValues)
            ;
    class_<DirichletMultinomialRobust>("DirichletMultinomialRobust")
            .constructor()
            .field("id", &DirichletMultinomialRobust::id)
            .method("SetLambdaValues", &DirichletMultinomialRobust::SetLambdaValues)
            ;
    class_<Multinomial>("Multinomial")
            .constructor()
            .field("id", &Multinomial::id)
            .method("SetLambdaValues", &Multinomial::SetLambdaValues)
            ;
    class_<MultinomialRobust>("MultinomialRobust")
            .constructor()
            .field("id", &MultinomialRobust::id)
            .method("SetLambdaValues", &MultinomialRobust::SetLambdaValues)
            ;
    class_<Fleet>("Fleet")
            .constructor()
            .field("id", &Fleet::id)
            .field("name", &Fleet::name)
            .field("catch_fraction_of_year", &Fleet::catch_fraction_of_year)
            .method("AddFishingMortality", &Fleet::AddFishingMortality)
            .method("AddSelectivity", &Fleet::AddSelectivity)
            .method("AddAgeCompData", &Fleet::AddAgeCompData)
            .method("AddIndexData", &Fleet::AddIndexData)
            .method("AddLengthCompData", &Fleet::AddLengthCompData)
            .method("SetAgeCompNllComponent", &Fleet::SetAgeCompNllComponent)
            .method("SetIndexNllComponent", &Fleet::SetIndexCompNllComponent)
            .method("SetLengthCompNllComponent", &Fleet::SetLengthCompNllComponent)
            ;

    class_<Survey>("Survey")
            .constructor()
            .field("id", &Survey::id)
            .field("name", &Survey::name)
            .field("survey_fraction_of_year", &Survey::survey_fraction_of_year)
            .field("q", &Survey::q)
            .method("AddSelectivity", &Survey::AddSelectivity)
            .method("AddAgeCompData", &Survey::AddAgeCompData)
            .method("AddIndexData", &Survey::AddIndexData)
            .method("AddLengthCompData", &Survey::AddLengthCompData)
            .method("SetAgeCompNllComponent", &Survey::SetAgeCompNllComponent)
            .method("SetIndexNllComponent", &Survey::SetIndexCompNllComponent)
            .method("SetLengthCompNllComponent", &Survey::SetLengthCompNllComponent)
            ;

    class_<IndexData>("IndexData")
            .constructor()
            .field("id", &IndexData::id)
            .field("values", &IndexData::data)
            .field("error", &IndexData::error)
            .field("sex", &IndexData::sex)
            .field("missing_value", &IndexData::missing_values)
            .field("is_abundance", &IndexData::is_abundance)
            ;

    class_<AgeCompData>("AgeCompData")
            .constructor()
            .field("id", &AgeCompData::id)
            .field("values", &AgeCompData::data)
            .field("sample_size", &AgeCompData::sample_size)
            .field("sex", &AgeCompData::sex)
            .field("missing_values", &AgeCompData::missing_values)
            ;
    class_<LengthCompData>("LengthCompData")
            .constructor()
            .field("id", &LengthCompData::id)
            .field("values", &LengthCompData::data)
            .field("sample_size", &LengthCompData::sample_size)
            .field("sex", &LengthCompData::sex)
            .field("missing_values", &LengthCompData::missing_values)
            ;

}

#endif /* R4MAS_HPP */

