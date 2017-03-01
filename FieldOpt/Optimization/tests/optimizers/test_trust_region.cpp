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
            Optimization::Optimizer *trust_region_search_ = new TrustRegionSearch(
                    settings_trust_region_search_max_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
            PolyModel poly_model = PolyModel(trust_region_search_->GetTentativeBestCase(), 1);
            //std::cout << "Created Polymodel" << std::endl;
            //poly_model.complete_points_abs();
            //std::cout << "Completed set of points" << std::endl;
            EXPECT_FALSE(poly_model.isModelReady());
            //poly_model.calculate_model_coeffs();
            //EXPECT_TRUE(poly_model.isModelReady());

        }


        // test Newton_Point
   TEST_F(TrustRegionSearchTest, Newton_Point) {
            test_case_2r_->set_objective_function_value(1000);
            Optimization::Optimizer *trust_region_search_ = new TrustRegionSearch(
                    settings_trust_region_search_max_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
            int dimension=2;
            int no_elemts_= (dimension+1)*(dimension+2)/2;
            Eigen::VectorXd points= Eigen::VectorXd::Random(dimension);
            Eigen::VectorXd coeffs_= Eigen::VectorXd::Random(no_elemts_);
            std::cout<<"Create polynomial:\n"<<dimension<<std::endl;
            std::cout<<"Dimension is:\n"<<dimension<<std::endl;
            std::cout<<"Point is:\n"<<points<<std::endl;
            std::cout<<"coefficent is:\n"<<coeffs_<<std::endl;
            Polynomial poly=Polynomial(dimension,coeffs_);
            poly.Hessian();
            std::cout<<"Find Quastion-Newton point:\n"<<std::endl;
            poly.Newton_Point(points);
            EXPECT_TRUE(true);


        }

        // sphere function, build model

        // test sphere
    TEST_F(TrustRegionSearchTest, Spherefunction) {
            test_case_2r_->set_objective_function_value(Sphere(test_case_2r_->GetRealVarVector()));
            Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_min_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
        Optimization::Case *tentative_best_0 = trust_region_search_->GetTentativeBestCase();
        for (int iter = 0; iter < 35; ++iter) {
            int No_of_case=iter+1;
            Optimization::Case *new_case = trust_region_search_->GetCaseForEvaluation();
            std::cout << "set objetive function value of case " << No_of_case<< std::endl;
            Eigen::VectorXd  point=new_case->GetRealVarVector();
            std::cout<<"point of this case is"<<point<<std::endl;
            new_case->set_objective_function_value(Sphere(new_case->GetRealVarVector()));
            std::cout << "the objective function valus is " << new_case->objective_function_value()<< std::endl;
            trust_region_search_->SubmitEvaluatedCase(new_case);
            trust_region_search_->optimizationStep();
        }


            Optimization::Case *tentative_best_final = trust_region_search_->GetTentativeBestCase();
            Eigen::VectorXd  point=tentative_best_final->GetRealVarVector();
            std::cout<<"point of final tentative_best_ case is"<<point<<std::endl;

            tentative_best_final->set_objective_function_value(Sphere(tentative_best_final->GetRealVarVector()));
            std::cout<<"objetive function value of final tentative_best_ case is"<<tentative_best_final->objective_function_value()<<std::endl;


            EXPECT_TRUE(tentative_best_final->objective_function_value() <= tentative_best_0->objective_function_value());
        }



        // Test matyas function for two variables
            TEST_F(TrustRegionSearchTest, Matyasfunction) {
                test_case_2r_->set_objective_function_value(Matyas(test_case_2r_->GetRealVarVector()));
                Optimization::Optimizer *trust_region_search_= new TrustRegionSearch(settings_trust_region_search_min_unconstr_, test_case_2r_, varcont_prod_bhp_, grid_5spot_);
                Optimization::Case *tentative_best_0 = trust_region_search_->GetTentativeBestCase();
                for (int iter = 0; iter < 11; ++iter) {
                    int No_of_case=iter+1;
                    Optimization::Case *new_case = trust_region_search_->GetCaseForEvaluation();
                    std::cout << "set objetive function value of case " << No_of_case<< std::endl;
                    Eigen::VectorXd  point=new_case->GetRealVarVector();
                    std::cout<<"point of this case is"<<point<<std::endl;
                    new_case->set_objective_function_value(Matyas(new_case->GetRealVarVector()));
                    std::cout << "the objective function valus is " << new_case->objective_function_value()<< std::endl;
                    trust_region_search_->SubmitEvaluatedCase(new_case);
                    trust_region_search_->optimizationStep();
                }


                Optimization::Case *tentative_best_final = trust_region_search_->GetTentativeBestCase();
                Eigen::VectorXd  point=tentative_best_final->GetRealVarVector();
                std::cout<<"point of final tentative_best_ case is"<<point<<std::endl;

                tentative_best_final->set_objective_function_value(Matyas(tentative_best_final->GetRealVarVector()));
                std::cout<<"objetive function value of final tentative_best_ case is"<<tentative_best_final->objective_function_value()<<std::endl;


                EXPECT_TRUE(tentative_best_final->objective_function_value() <= tentative_best_0->objective_function_value());
            }
}

