#include <gtest/gtest.h>
#include "Optimization/tests/test_resource_optimizer.h"
#include "../Runner/bookkeeper.h"
#include "Optimization/optimizers/compass_search.h"

namespace {

    class BookkeeperTest : public ::testing::Test, public TestResources::TestResourceOptimizer {
    protected:
        BookkeeperTest() {
            test_case_2r_->set_objective_function_value(100.0);
            compass_search_ = new ::Optimization::Optimizers::CompassSearch(settings_compass_search_max_unconstr_,
                                                                            test_case_2r_, model_->variables(),
                                                                            model_->grid());
            bookkeeper_ = new Runner::Bookkeeper(settings_full_, compass_search_->case_handler());
            c1 = compass_search_->GetCaseForEvaluation();
            c2 = compass_search_->GetCaseForEvaluation();
            c3 = compass_search_->GetCaseForEvaluation();
        }
        virtual ~BookkeeperTest() {}
        virtual void SetUp() {}

        Runner::Bookkeeper *bookkeeper_;
        Optimization::Optimizer *compass_search_;
        Optimization::Case *c1;
        Optimization::Case *c2;
        Optimization::Case *c3;
    };

    TEST_F(BookkeeperTest, Bookkeeping) {
        EXPECT_FALSE(base_case_->Equals(c1));
        EXPECT_FALSE(base_case_->Equals(c2));
        EXPECT_FALSE(base_case_->Equals(c3));
        EXPECT_TRUE(test_case_2r_->Equals(test_case_2r_));
        EXPECT_TRUE(bookkeeper_->IsEvaluated(test_case_2r_));
        EXPECT_FALSE(bookkeeper_->IsEvaluated(c1));
        EXPECT_FALSE(bookkeeper_->IsEvaluated(c2));
        EXPECT_FALSE(bookkeeper_->IsEvaluated(c3));

        c1->set_objective_function_value(100);
        compass_search_->SubmitEvaluatedCase(c1);
        EXPECT_TRUE(bookkeeper_->IsEvaluated(c1));
    }

}
