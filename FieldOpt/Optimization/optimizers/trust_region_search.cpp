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
           // std::cout << "Radius before scaling is  " << radius_ << std::endl;
            std::cout << " Factor k is  " << k<< std::endl;
            radius_ = k*radius_;
            polymodel_.setRadius(radius_);
            std::cout << "Radius after scaling is " << radius_ << std::endl;
        }

        //tentative_best_case_=base_case(Optimizer(settings,base_case,variable,grid)
        // case_handler_ = new CaseHandler(tentative_best_case_);

        Optimizer::TerminationCondition TrustRegionSearch::IsFinished()
        {  std::cout << "EvaluatedCases().size is\n  " <<case_handler_->EvaluatedCases().size()<<std::endl;
         /*   std::cout << "Radius is\n" <<radius_<<std::endl;*/
            if (case_handler_->EvaluatedCases().size() >= max_evaluations_)
                return MAX_EVALS_REACHED;
            else if (radius_ < minimum_radius_)
                return MINIMUM_STEP_LENGTH_REACHED;
            else return NOT_FINISHED;

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

            need_optimization_step=true;
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
        }

        void TrustRegionSearch::UpdateModel() {
            if (!need_optimization_step)
            {
                if(case_handler_->QueuedCases().size()==0)
                {
                    std::cout <<"The new Base case has been evaluated after optimization step, we need to check the model accuracy" << std::endl;
                    double f_xk=currentBaseCase->objective_function_value();
                    double f_xk1=newBaseCase->objective_function_value();
                    double m_xk=polymodel_.obejctive_function_value_model(Current_CenterPoint);
                    double m_xk1=polymodel_.obejctive_function_value_model(New_CenterPoint);
                    double ared=f_xk-f_xk1;
                    double pred=m_xk-m_xk1;
                    double rho=ared/pred;
                    std::cout << "The objective function value of current base case from model is\n " <<m_xk << std::endl;
                    std::cout << "The objective function value of current new base case from model is\n " <<m_xk1<< std::endl;
                    std::cout <<"Actual objective function value of the current Base Case is \n" <<f_xk<<std::endl;
                    std::cout <<"Actual objective function value of the new Base Case is \n" <<f_xk1<<std::endl;
                    std::cout <<"Actual reduction is  " << ared<<std::endl;
                    std::cout <<"Model reduction is " << pred<<std::endl;
                    Optimization::Case *BaseCase;
                    Eigen::VectorXd  CenterPoint;

                   /*     std::cout <<"f(xk) is " <<f_xk<<std::endl;
                   std::cout <<"f(xk+1) is " <<f_xk1<<std::endl;
                   std::cout <<"m(xk) is " <<m_xk<<std::endl;
                   std::cout <<"m(xk1) is " <<m_xk1<<std::endl;

                   */
                    if(rho>=0.75)
                    {    std::cout <<"ratio is " << rho<<", ratio >= 0.75"<<std::endl;
                         std::cout <<"Very successful iteration and candidate solution is accepted, radius should be increased for next iteration " <<std::endl;
                        scaleRadius(2);
                         BaseCase=newBaseCase;
                         CenterPoint=New_CenterPoint;
                    }
                    else if(rho>=0.25)
                    {   std::cout <<"ratio is " << rho<<", 0.25<=ratio<=0.75"<<std::endl;
                        std::cout <<"Successful iteration and candidate solution is accpeted,keep same radius for next iteration " << std::endl;
                        BaseCase=newBaseCase;
                        CenterPoint=New_CenterPoint;
                    }

                    else
                    {    std::cout <<"ratio is " << rho<<", ratio<=0.25"<<std::endl;
                        std::cout <<"unsuccessful iteration and candidate solution is not accpeted, radius should be reduced for next iteration " << std::endl;
                        scaleRadius(0.5);
                        BaseCase=currentBaseCase;
                        CenterPoint=Current_CenterPoint;

                    }


                    polymodel_.addBaseCase(BaseCase); //needs_set of points=true, is model_complete=false;
                    tentative_best_case_ =BaseCase;
                    std::cout <<"CenterPoint for next itertion is \n" << CenterPoint<<std::endl;

                }
                else
                {
                    std::cout <<" Optimization step is finished, the new base case needs to be evaluated first." << std::endl;
                }

            }
        }

        void TrustRegionSearch::optimizationStep() //
        {
            if (case_handler_->QueuedCases().size()==0&&need_optimization_step) {
                std::cout << "All cases have been evaluated, we make optimization step now " << std::endl;

                polymodel_.calculate_model_coeffs();
                Eigen::VectorXd optimizationstep;

                if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
                    optimizationstep = -polymodel_.optimizationStep_CP();
                   // optimizationstep = -polymodel_.optimizationStep_SDL();

                    std::cout << "Optimizer Mode is Maximize, the optimization step should be   " << optimizationstep
                              << std::endl;
                } else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
                    optimizationstep = polymodel_.optimizationStep_CP();
                    //optimizationstep = polymodel_.optimizationStep_SDL();
                    std::cout << "Optimizer Mode is minimize, The optimization step shoule be  " << optimizationstep
                              << std::endl;
                }
                std::cout << "The current Center Point is \n" << polymodel_.get_centerpoint() << std::endl;
                 currentBaseCase=tentative_best_case_;
                 Current_CenterPoint=polymodel_.get_centerpoint();
                New_CenterPoint = optimizationstep + polymodel_.get_centerpoint();
                std::cout << "The New Center Point after optimization step is\n" << New_CenterPoint << std::endl;
                // Use the best case in current subregion as a new Base Case for next iteration. Creat this Case from its Point
                // tentative_best_case=base_case (center point)
                newBaseCase = polymodel_.CaseFromPoint(New_CenterPoint, tentative_best_case_);
                case_handler_->AddNewCase(newBaseCase);

                need_optimization_step=false;

            }
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

        void TrustRegionSearch::SubmitEvaluatedCase(Case *c) {
            
            case_handler_->UpdateCaseObjectiveFunctionValue(c->id(), c->objective_function_value());
            case_handler_->SetCaseEvaluated(c->id());
            handleEvaluatedCase(case_handler_->GetCase(c->id()));
            optimizationStep();
            UpdateModel();
            
        }

    }}

