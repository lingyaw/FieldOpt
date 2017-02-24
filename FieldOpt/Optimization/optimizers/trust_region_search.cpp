#include "trust_region_search.h"

namespace Optimization {
    namespace Optimizers {

        TrustRegionSearch::TrustRegionSearch(::Settings::Optimizer *settings, Case *base_case,
                                     Model::Properties::VariablePropertyContainer *variables,
                                     Reservoir::Grid::Grid *grid)
                : Optimizer(settings, base_case, variables, grid)
        {
            radius_ = settings->parameters().initial_step_length;
            minimum_radius_ = settings->parameters().minimum_step_length;
        }

        void TrustRegionSearch::scaleRadius(double k)
        {
            radius_ = k*radius_;
            polymodel_.setRadius(k);
        }

        //tentative_best_case_=base_case(Optimizer(settings,base_case,variable,grid)
        // case_handler_ = new CaseHandler(tentative_best_case_);

        Optimizer::TerminationCondition TrustRegionSearch::IsFinished()
        {
            if (case_handler_->EvaluatedCases().size() >= max_evaluations_)
                return MAX_EVALS_REACHED;
            else if (radius_ < minimum_radius_)
                return MINIMUM_STEP_LENGTH_REACHED;
            else return NOT_FINISHED; // The value of not finished is 0, which evaluates to false.
        }

        void TrustRegionSearch::iterate() {
            std::cout << "run iterate() and create interpolation points " << std::endl;
            std::cout << "this is iteration number " << iteration_ << std::endl;
            /* Everytime we update model we must first have a PolyModel object,
             * then we must complete the set of points in the model, then the
             * objective values of all cases must be calculated, and lastly we
             * can return the PolyModel coefficients. This functions checks these
             * steps, starting with the first one. We start at the first incomplete
             * step and perform all the steps that come after it.
             */

            // At the first iteration we initialze the PolyModel with base_case
            if (iteration_ == 0) {
                std:: cout << "Initializing polymodel and completing points" << std::endl;
                initializeModel();
                std::cout <<"Is model completed  after initializeModel? "<<polymodel_.isModelReady() << std::endl;


            }

            /* Either the PolyModel is ready to give us a derivative, or
             * we need to update points or get cases evaluated
             */
            else if(!polymodel_.isModelReady()) {
                std:: cout << "Model is not ready, Completing model first" << std::endl;
                completeModel();

            }

                // if all cases have been evaluated(set objective value), calculate moodel coefficient, mode is ready

                // if (-----){
                //
                // polymodel_.calculate_model_coeffs();
                //
                // }


            else {

                std::cout << "Is model completed (else)? " << polymodel_.isModelReady() << std::endl;
                std::cout << "Model should be ready, we make some opt step " << std::endl;
                //TODO Here we can use current_model to do an optimization step, which
                optimizationStep();
                std::cout << "end of opt step part" << std::endl;
            }

            case_handler_->ClearRecentlyEvaluatedCases();
            iteration_++;
        }


        void TrustRegionSearch::initializeModel() {
            // Create polymodel with radius and center and complete set of points.
            polymodel_ = PolyModel(GetTentativeBestCase(), radius_); // is_model_complete=false; need_set_of points=true;
            completeModel();
        }

        void TrustRegionSearch::completeModel() {
            polymodel_.complete_points(); // needs_set of points=false after completing points
            // Add cases to case_handler and clear CasesNotEval queue
            case_handler_->AddNewCases(polymodel_.get_cases_not_eval());
            polymodel_.ClearCasesNotEval(); // needs_evals=false
            polymodel_.set_model_complete(); //is_model_complete=true
            polymodel_.set_evaluations_complete();
        }

        void TrustRegionSearch::optimizationStep() // we dont need Perturbation, gradient descent og dog-leg method
        {
            polymodel_.calculate_model_coeffs();
            Eigen::VectorXd optimizationstep;

            /* optimizationstep from polymodel for minimum value.
             * Optimization step need to be for maximum value
             * check model_( maximum or minimum)
            */

            if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
              optimizationstep=-polymodel_.optimizationStep_NG();
                std::cout << "Optimizer Mode is Maximize, the optimization step should be   " << optimizationstep<< std::endl;
            }
            else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
                optimizationstep=polymodel_.optimizationStep_NG();
                std::cout << "Optimizer Mode is minimize, The optimization step shoule be  " << optimizationstep<< std::endl;
            }
            // Find point of newCenterPoint which is the best point in current subregion
            // newCenterpoint=current centerpoint+optimization step
            Eigen::VectorXd  NewCenterPoint=optimizationstep+polymodel_.get_centerpoint();
            std::cout << "The current Center Point is   " << polymodel_.get_centerpoint()<< std::endl;
            std::cout << "The New Center Point is   " << NewCenterPoint<< std::endl;

            // Use the best case in current subregion as a new Base Case for next iteration. Creat this Case from its Point
            // tentative_best_case=base_case (center point)
            Case *newBaseCase=polymodel_.CaseFromPoint(NewCenterPoint,tentative_best_case_);
            polymodel_.addCenterPoint(NewCenterPoint); //needs_set of points=true, is model_complete=false;
            std::cout <<"Is model completed after add CenterPoint? "<<polymodel_.isModelReady() << std::endl;
            scaleRadius(0.5);
            tentative_best_case_=newBaseCase;
            //Eigen::VectorXd tentative_best_case=tentative_best_case_->GetRealVarVector();
           // std::cout << "tentative best case is (New Center Point)  " << tentative_best_case<< std::endl;



            //polymodel_.optimizationStep();
            //Case* a =polymodel_.find_NewBaseCase();
            // applyNew.. sets best case so far to new best.
            // let's just go with the flow for now
            //tentative_best_case=base_case;
            //QList<Case *> perturbations = QList<Case *>();
            //for (QUuid id : tentative_best_case_->integer_variables().keys())
             //   perturbations.append(tentative_best_case_->Perturb(id, Case::SIGN::PLUS, radius_/3));
           // Case* c = perturbations.at(0);

            //case_handler_->AddNewCase(c); //AddNewCase Add a new non-evaluated case to the evaluation_queue_.



            //? tentative_best_case=c(new center point)
           // applyNewTentativeBestCase(); // NewTentativeBestCase from RecentlyEvaluatedCases
            /* Add new center point, this also sets
             * polymodel_isModelReady to false.
             */





        }

        QString TrustRegionSearch::GetStatusStringHeader() const
        {
            return QString("%1,%2,%3,%4,%5,%6,%7")
                    .arg("Iteration")
                    .arg("EvaluatedCases")
                    .arg("QueuedCases")
                    .arg("RecentlyEvaluatedCases")
                    .arg("TentativeBestCaseID")
                    .arg("TentativeBestCaseOFValue")
                    .arg("StepLength");
        }

        QString TrustRegionSearch::GetStatusString() const
        {
            return QString("%1,%2,%3,%4,%5,%6,%7")
                    .arg(iteration_)
                    .arg(nr_evaluated_cases())
                    .arg(nr_queued_cases())
                    .arg(nr_recently_evaluated_cases())
                    .arg(tentative_best_case_->id().toString())
                    .arg(tentative_best_case_->objective_function_value())
                    .arg(radius_);
        }

    }}

