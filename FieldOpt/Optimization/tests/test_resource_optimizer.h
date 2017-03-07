//
// Created by einar on 6/2/16.
//

#ifndef FIELDOPT_TEST_RESOURCE_OPTIMIZER_H
#define FIELDOPT_TEST_RESOURCE_OPTIMIZER_H

#include "Model/tests/test_resource_model.h"
#include "Optimization/case.h"
#include "test_resource_cases.h"

namespace TestResources {
    class TestResourceOptimizer : public TestResourceModel, public TestResourceCases {
    protected:
        TestResourceOptimizer() {
            base_case_ = new Optimization::Case(model_->variables()->GetBinaryVariableValues(),
                                                model_->variables()->GetDiscreteVariableValues(),
                                                model_->variables()->GetContinousVariableValues());
            base_case_->set_objective_function_value(1000.0);

            settings_compass_search_min_unconstr_ = new Settings::Optimizer(get_json_settings_compass_search_minimize_);
            settings_compass_search_max_unconstr_ = new Settings::Optimizer(get_json_settings_compass_search_maximize_);
            settings_apps_min_unconstr_ = new Settings::Optimizer(get_json_settings_apps_minimize_);
            settings_apps_max_unconstr_ = new Settings::Optimizer(get_json_settings_apps_maximize_);
            settings_trust_region_search_min_unconstr_ = new Settings::Optimizer(get_json_settings_trust_region_search_minimize_);
            settings_trust_region_search_max_unconstr_ = new Settings::Optimizer(get_json_settings_trust_region_search_maximize_);

        }

        Optimization::Case *base_case_;
        Settings::Optimizer *settings_compass_search_min_unconstr_;
        Settings::Optimizer *settings_compass_search_max_unconstr_;
        Settings::Optimizer *settings_apps_min_unconstr_;
        Settings::Optimizer *settings_apps_max_unconstr_;
        Settings::Optimizer *settings_trust_region_search_min_unconstr_;
        Settings::Optimizer *settings_trust_region_search_max_unconstr_;

    private:
        QJsonObject obj_fun_ {
                {"Type", "WeightedSum"},
                {"WeightedSumComponents", QJsonArray{
                        QJsonObject{
                                {"Coefficient", 1.0}, {"Property", "CumulativeOilProduction"}, {"TimeStep", -1}, {"IsWellProp", false}
                        },
                        QJsonObject{
                                {"Coefficient", 0.0}, {"Property", "CumulativeWaterProduction"}, {"TimeStep", -1}, {"IsWellProp", false}
                        }
                }}
        };
        QJsonObject get_json_settings_compass_search_minimize_ {
                {"Type", "Compass"},
                {"Mode", "Minimize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 100},
                        {"InitialStepLength", 0.25},
                        {"MinimumStepLength", 0.01}
                }},
                {"Objective", obj_fun_}
        };

        QJsonObject get_json_settings_compass_search_maximize_ {
                {"Type", "Compass"},
                {"Mode", "Maximize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 100},
                        {"InitialStepLength", 8},
                        {"MinimumStepLength", 1}
                }},
                {"Objective", obj_fun_}
        };

        QJsonObject get_json_settings_apps_minimize_ {
                {"Type", "APPS"},
                {"Mode", "Minimize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 500},
                        {"InitialStepLength", 0.64},
                        {"MinimumStepLength", 0.005}
                }},
                {"Objective", obj_fun_}
        };

        QJsonObject get_json_settings_apps_maximize_ {
                {"Type", "APPS"},
                {"Mode", "Maximize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 500},
                        {"InitialStepLength", 8},
                        {"MinimumStepLength", 1}
                }},
                {"Objective", obj_fun_}
        };
        QJsonObject get_json_settings_trust_region_search_minimize_ {
                {"Type", "Trustregion"},
                {"Mode", "Minimize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 100},
                        {"InitialStepLength", 1},
                        {"MinimumStepLength", 0.005}
                }},
                {"Objective", obj_fun_}
        };
        QJsonObject get_json_settings_trust_region_search_maximize_ {
                {"Type", "Trustregion"},
                {"Mode", "Maximize"},
                {"Parameters", QJsonObject{
                        {"MaxEvaluations", 100},
                        {"InitialStepLength", 8},
                        {"MinimumStepLength", 0.005}
                }},
                {"Objective", obj_fun_}
        };


    };
}

#endif //FIELDOPT_TEST_RESOURCE_OPTIMIZER_H
