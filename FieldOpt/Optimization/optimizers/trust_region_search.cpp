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
            std::cout << "Radius before scaling is  " << radius_ << std::endl;
            std::cout << " Factor k is  " << k<< std::endl;
            radius_ = k*radius_;
            polymodel_.setRadius(radius_);
            std::cout << "Radius after scaling is " << radius_ << std::endl;
        }

        //tentative_best_case_=base_case(Optimizer(settings,base_case,variable,grid)
        // case_handler_ = new CaseHandler(tentative_best_case_);

        Optimizer::TerminationCondition TrustRegionSearch::IsFinished()
        {/*   std::cout<<"EvauatedCases().size()"<<case_handler_->EvaluatedCases().size()<<std::endl;
            std::cout<<"queuedCases().size()"<<case_handler_->QueuedCases().size()<<std::endl;
            std::cout<<"grad.norm"<<grad_norm<<std::endl;*/
            if (case_handler_->EvaluatedCases().size() >= max_evaluations_)
                return MAX_EVALS_REACHED;
            else if (radius_ < minimum_radius_)
                return MINIMUM_STEP_LENGTH_REACHED;
            else return NOT_FINISHED; // The value of not finished is 0, which evaluates to false.
            /*
            if (case_handler_->EvaluatedCases().size() >= max_evaluations_)
                return MAX_EVALS_REACHED;
            else  (radius_ < minimum_radius_)
                return MINIMUM_STEP_LENGTH_REACHED;
            //else if (iteration_ == 0) return NOT_FINISHED;
            //else if (case_handler_->QueuedCases().size() > 0 ) return NOT_FINISHED;
            //else if (grad_norm>=epsilon) return NOT_FINISHED;
           // else return MAX_EVALS_REACHED;*/


        }



        void TrustRegionSearch::iterate() {
            std::cout << "Run iterate() and create interpolation points " << std::endl;
            std::cout << "This is iteration:\n " << iteration_ << std::endl;
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


            }

            /* Either the PolyModel is ready to give us a derivative, or
             * we need to update points or get cases evaluated
             */
            else  {
                std:: cout << "Model is not ready, Completing model first" << std::endl;
                completeModel();

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
            //polymodel_.set_model_complete(); //is_model_complete=true
            polymodel_.set_evaluations_complete();
            need_optimization_step= true;
        }
/*
        void TrustRegionSearch::optimizationStep() // we dont need Perturbation, gradient descent og dog-leg method
        {
            // all cases must been evaluted first befor optimization step
            if (case_handler_->QueuedCases().size()==0) {
                std::cout << "All cases have been evaluated, we make optimization step now " << std::endl;

                polymodel_.calculate_model_coeffs();
                grad_norm=polymodel_.Gradient().norm();
                Eigen::VectorXd optimizationstep;

                *//* optimizationstep from polymodel for minimum value.
                 * Optimization step need to be for maximum value
                 * check model_( maximum or minimum)
                *//*

                if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
                   // optimizationstep = -polymodel_.optimizationStep_CP();
                    optimizationstep = -polymodel_.optimizationStep_SDL();

                    std::cout << "Optimizer Mode is Maximize, the optimization step should be   " << optimizationstep
                              << std::endl;
                } else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
                   // optimizationstep = polymodel_.optimizationStep_CP();
                    optimizationstep = polymodel_.optimizationStep_SDL();
                    std::cout << "Optimizer Mode is minimize, The optimization step shoule be  " << optimizationstep
                              << std::endl;
                }
                // Find point of newCenterPoint which is the best point in current subregion
                // newCenterpoint=current centerpoint+optimization step
                std::cout << "The current Center Point is \n" << polymodel_.get_centerpoint() << std::endl;
                std::cout << "The objective function value of current center is\n " << polymodel_.obejctive_function_value_model(polymodel_.get_centerpoint()) << std::endl;
               // currentBaseCase=tentative_best_case_;
               // Current_CenterPoint=polymodel_.get_centerpoint();
                New_CenterPoint = optimizationstep + polymodel_.get_centerpoint();
                std::cout << "The New Center Point after optimization step is\n" << New_CenterPoint << std::endl;

                // Use the best case in current subregion as a new Base Case for next iteration. Creat this Case from its Point
                // tentative_best_case=base_case (center point)
                newBaseCase = polymodel_.CaseFromPoint(New_CenterPoint, tentative_best_case_);


                polymodel_.addCenterPoint(New_CenterPoint); //needs_set of points=true, is model_complete=false;
                //scaleRadius(0.5);
                //handleEvaluatedCase(newBaseCase);
                tentative_best_case_ = newBaseCase;
                //std::cout << "points of tentative_best_case is   " << tentative_best_case_->GetRealVarVector() << std::endl;

                //calulate the objective function value of new center point (based on polymodel)
                objective_value = polymodel_.obejctive_function_value_model(tentative_best_case_->GetRealVarVector());

               std::cout << "The objective function value of new center point based on polymodel is "
                          << objective_value << std::endl;
                std::cout << "|| g|| at current center point is: " << grad_norm<<std::endl;
                if(grad_norm<=epsilon){
                    std::cout << "|| g|| <="<<epsilon<<". This is the final iteration " <<std::endl;
                }
                else
                {
                    std::cout << "|| g|| >"<<epsilon<<". We should contiue to next iteration"<<std::endl;
                }

            }


        }*/


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

        void TrustRegionSearch::SubmitEvaluatedCase(Case *c) {
            
            case_handler_->UpdateCaseObjectiveFunctionValue(c->id(), c->objective_function_value());
            case_handler_->SetCaseEvaluated(c->id());
            handleEvaluatedCase(case_handler_->GetCase(c->id()));

            // all cases must been evaluted first before optimization step
            if (case_handler_->QueuedCases().size()==0) {
                std::cout << "All cases have been evaluated, we make optimization step now " << std::endl;

                polymodel_.calculate_model_coeffs();
               // grad_norm=polymodel_.Gradient().norm();
                Eigen::VectorXd optimizationstep;

                /* optimizationstep from polymodel for minimum value.
                 * Optimization step need to be for maximum value
                 * check model_( maximum or minimum)
                */

                if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
                   //optimizationstep = -polymodel_.optimizationStep_CP();
                   optimizationstep = -polymodel_.optimizationStep_SDL();

                    std::cout << "Optimizer Mode is Maximize, the optimization step should be   " << optimizationstep
                              << std::endl;
                } else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
                    //optimizationstep = polymodel_.optimizationStep_CP();
                    optimizationstep = polymodel_.optimizationStep_SDL();
                    std::cout << "Optimizer Mode is minimize, The optimization step shoule be  " << optimizationstep
                              << std::endl;
                }
                // Find point of newCenterPoint which is the best point in current subregion
                // newCenterpoint=current centerpoint+optimization step
                std::cout << "The current Center Point is \n" << polymodel_.get_centerpoint() << std::endl;
                std::cout << "The objective function value of current center is\n " << polymodel_.obejctive_function_value_model(polymodel_.get_centerpoint()) << std::endl;
                // currentBaseCase=tentative_best_case_;
                // Current_CenterPoint=polymodel_.get_centerpoint();
                New_CenterPoint = optimizationstep + polymodel_.get_centerpoint();
                std::cout << "The New Center Point after optimization step is\n" << New_CenterPoint << std::endl;

                // Use the best case in current subregion as a new Base Case for next iteration. Creat this Case from its Point
                // tentative_best_case=base_case (center point)
                newBaseCase = polymodel_.CaseFromPoint(New_CenterPoint, tentative_best_case_);


                polymodel_.addCenterPoint(New_CenterPoint); //needs_set of points=true, is model_complete=false;
                //scaleRadius(0.5);
                //handleEvaluatedCase(newBaseCase);
                tentative_best_case_ = newBaseCase;
                //std::cout << "points of tentative_best_case is   " << tentative_best_case_->GetRealVarVector() << std::endl;

                //calulate the objective function value of new center point (based on polymodel)
                objective_value = polymodel_.obejctive_function_value_model(tentative_best_case_->GetRealVarVector());

                std::cout << "The objective function value of new center point based on polymodel is "
                          << objective_value << std::endl;

              /* // std::cout << "|| g|| at current center point is: " << grad_norm<<std::endl;
                grad_norm=polymodel_.Gradient().norm();
                std::cout << "|| g|| at new center point is: " << grad_norm<<std::endl;

                if(grad_norm<=epsilon){
                    std::cout << "|| g|| <="<<epsilon<<". This is the final iteration " <<std::endl;
                }
                else
                {
                    std::cout << "|| g|| >"<<epsilon<<". We should contiue to next iteration"<<std::endl;
                }*/

            }
            
        }

    }}

