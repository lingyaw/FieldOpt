#ifndef FIELDOPT_TEST_RESOURCE_CASES_H
#define FIELDOPT_TEST_RESOURCE_CASES_H

#include "Optimization/case.h"
#include "Model/properties/continous_property.h"
#include "Model/tests/test_resource_variable_property_container.h"
#include <QHash>
#include <QUuid>
#include <QList>

namespace TestResources {
    class TestResourceCases : public TestResources::TestResourceVariablePropertyContainer {
    public:
        TestResourceCases() {
            test_case_1_3i_ = new Optimization::Case(QHash<QUuid, bool>(), integer_variables_3d_, QHash<QUuid, double>());
            test_case_2_3r_ = new Optimization::Case(QHash<QUuid, bool>(), QHash<QUuid, int>(), real_variables_3d_);
            test_case_2_3r_->set_objective_function_value(100.0);
            test_case_3_4b3i3r_ = new Optimization::Case(binary_variables_4d_, integer_variables_3d_, real_variables_3d_);
            test_case_3_4b3i3r_->set_objective_function_value(-50.0);
            test_case_4_4b3i3r = new Optimization::Case(binary_variables_4d_, integer_variables_3d_, real_variables_3d_);
            test_case_4_4b3i3r->set_objective_function_value(-50.0);
            trivial_cases_ << test_case_1_3i_ << test_case_2_3r_ << test_case_3_4b3i3r_ << test_case_4_4b3i3r;

            test_case_spline_ = new Optimization::Case(varcont_prod_spline_->GetBinaryVariableValues(),
                                                       varcont_prod_spline_->GetDiscreteVariableValues(),
                                                       varcont_prod_spline_->GetContinousVariableValues());

            test_case_2r_ = new Optimization::Case(QHash<QUuid, bool>(), QHash<QUuid, int>(), real_variables_2d_);
            test_case_2_4r_ = new Optimization::Case(QHash<QUuid, bool>(), QHash<QUuid, int>(), real_variables_4d_);
            test_case_2_19r_ = new Optimization::Case(QHash<QUuid, bool>(), QHash<QUuid, int>(), real_variables_19d_);
        }

        QList<Optimization::Case *> trivial_cases_;

        /* Case 1:
         * Only integer variables.
         * OBjective function not evaluated.
         */
        Optimization::Case *test_case_1_3i_;

        /* Case 2:
         * Only real variables.
         * Objective function evaluated.
         */
        Optimization::Case *test_case_2_3r_;

        /* Case 3:
         * All variable types.
         * Objective function evaluated.
         */
        Optimization::Case *test_case_3_4b3i3r_;

        /* Case 4:
         * Identical to case 3.
         */
        Optimization::Case *test_case_4_4b3i3r;

        /* Case 5:
         * Spline well case
         */
        Optimization::Case *test_case_spline_;

        /* Case:
         * Two real variables.
         */
        Optimization::Case *test_case_2r_;
        /* Case:
        * Four real variables.
        */
        Optimization::Case *test_case_2_4r_;
        Optimization::Case *test_case_2_19r_;


    private:
        const QHash<QUuid, bool> binary_variables_4d_{
                {QUuid::createUuid(), true},
                {QUuid::createUuid(), true},
                {QUuid::createUuid(), false},
                {QUuid::createUuid(), false}
        };
        const QHash<QUuid, int> integer_variables_3d_{
                {QUuid::createUuid(), 1},
                {QUuid::createUuid(), 2},
                {QUuid::createUuid(), 5}
        };
        const QHash<QUuid, double> real_variables_19d_{
                {QUuid::createUuid(), 1.0},
                {QUuid::createUuid(), 1.0},
                {QUuid::createUuid(), 100.0},
                {QUuid::createUuid(), 500.0},
                {QUuid::createUuid(), 100.0},
                {QUuid::createUuid(), 36.0},
                {QUuid::createUuid(), 1712.0},
                {QUuid::createUuid(), 12.0},
                {QUuid::createUuid(), 2000.0},
                {QUuid::createUuid(), 12.0},
                {QUuid::createUuid(), 2000.0},
                {QUuid::createUuid(), 1.0},
                {QUuid::createUuid(), 12.0},
                {QUuid::createUuid(), 500.0},
                {QUuid::createUuid(), 1712.0},
                {QUuid::createUuid(), 1712.0},
                {QUuid::createUuid(), 2000.0},
                {QUuid::createUuid(), 2000.0},
                {QUuid::createUuid(), 1712.0}

        };
        const QHash<QUuid, double> real_variables_4d_{
                {QUuid::createUuid(), 1.0},
                {QUuid::createUuid(), 4.0},
                {QUuid::createUuid(), 2.5},
                {QUuid::createUuid(), 5.5}
        };
        const QHash<QUuid, double> real_variables_3d_{
                {QUuid::createUuid(), 1.0},
                {QUuid::createUuid(), 4.0},
                {QUuid::createUuid(), 2.5}
        };

        const QHash<QUuid, double> real_variables_2d_{
                {QUuid::createUuid(), 4},
                {QUuid::createUuid(), 5}
        };

    };

}
#endif //FIELDOPT_TEST_RESOURCE_CASES_H
