#include <gtest/gtest.h>
#include "Optimization/optimizers/trust_region_search.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
    namespace {


    class TrustRegionSearchTest : public ::testing::Test,
                                  public TestResources::TestResourceOptimizer,
                                  public TestResources::TestResourceGrids {
    protected:
        TrustRegionSearchTest() {
            base_=base_case_;
            //trust_region_search_ = new ::Optimization::Optimizers::TrustRegionSearch(settings_optimizer_, base_case_, model_->variables(), grid_5spot_);
        }
        virtual ~TrustRegionSearchTest() {}
        virtual void SetUp() {}
        Optimization::Case *base_;
        //Optimization::Optimizer *trust_region_search_;
    };

    TEST_F(TrustRegionSearchTest, Constructor) {
            test_case_2r_->set_objective_function_value(1000);
            Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_min_unconstr_,
                                                                                 test_case_2r_,
                                                                                 varcont_prod_bhp_,
                                                                                 grid_5spot_);

        EXPECT_FLOAT_EQ(1000.0, trust_region_search_->GetTentativeBestCase()->objective_function_value());
    }

    TEST_F(TrustRegionSearchTest, BogoTest) {
        EXPECT_TRUE(true);
    }

    TEST_F(TrustRegionSearchTest, PointCaseConversionTest) {
            test_case_2r_->set_objective_function_value(1000);
            Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_max_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
        Optimization::Case* a = trust_region_search_->GetTentativeBestCase();
        Eigen::VectorXd vec = a->GetRealVarVector();
        Optimization::Case* b = PolyModel::CaseFromPoint(a->GetRealVarVector(), a);
        EXPECT_TRUE(a->Equals(b));
    }

    TEST_F(TrustRegionSearchTest, CompleteSetOfPointsTest) {
            test_case_2r_->set_objective_function_value(1000);
            Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_max_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
        PolyModel poly_model = PolyModel(trust_region_search_->GetTentativeBestCase(), 1);
        std::cout << "Created Polymodel" << std::endl;
        poly_model.complete_points();
        std::cout << "Completed set of points" << std::endl;
        EXPECT_FALSE(poly_model.isModelReady());
        //poly_model.calculate_model_coeffs();
        //EXPECT_TRUE(poly_model.isModelReady());
    }

        // sphere function, build model

        // test sphere
    TEST_F(TrustRegionSearchTest, OneIterationTest) {
            test_case_2r_->set_objective_function_value(Sphere(test_case_2r_->GetRealVarVector()));
            Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_min_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
        Optimization::Case *tentative_best_0 = trust_region_search_->GetTentativeBestCase();
        for (int iter = 0; iter <6; ++iter) {
            int No_of_case=iter+1;
            Optimization::Case *new_case = trust_region_search_->GetCaseForEvaluation();
            std::cout << "set objetive function value of case " << No_of_case<< std::endl;
            new_case->set_objective_function_value(Sphere(new_case->GetRealVarVector()));
            std::cout << "the objective function valus is " << new_case->objective_function_value()<< std::endl;
            trust_region_search_->SubmitEvaluatedCase(new_case);
        }
           //trust_region_search_->optimizationStep();

            Optimization::Case *tentative_best_1 = trust_region_search_->GetTentativeBestCase();
            // set objective function value of tentative_best_1
             //double objective_value=trust_region_search_->objective_value_model();
          ///  tentative_best_1->set_objective_function_value(objective_value);


            EXPECT_TRUE(tentative_best_1->objective_function_value() == tentative_best_0->objective_function_value());
    }
}

