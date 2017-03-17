//
// Created by cutie on 07.10.16.
//

#include "poly_model.h"

PolyModel::PolyModel(Optimization::Case* initial_case, double radius) {
    cases_.append(initial_case);
    points_.append(initial_case->GetRealVarVector());
    center_ = points_.at(0);
    std::cout << "Center point when creating model is:\n " << center_ << std::endl;
    radius_ = radius;
    std::cout << "Radius when creating model is: \n" << radius_ << std::endl;
    std::cout << "Objective value of current center point is: \n " << initial_case->objective_function_value()<< std::endl;

    dimension_ = center_.size();
    QList<Polynomial> basis;
    for (int i = 0; i < (dimension_+1)*(dimension_+2)/2; ++i) {
        Eigen::VectorXd temp_vec = Eigen::VectorXd::Zero((dimension_+1)*(dimension_+2)/2);
        temp_vec(i) = 1;
        Polynomial temp_poly = Polynomial(dimension_, temp_vec);
        basis.append(temp_poly);
    }
    basis_ = basis;

    needs_set_of_points_ = true;

  /*  // Check if objective value of base case has been computed
    if(initial_case->objective_function_value() == std::numeric_limits<double>::max()) {
        needs_evals_ = true;
        cases_not_eval_.append(initial_case);
    }
    else{needs_evals_ = false;}*/
    needs_evals_ = false;
    is_model_complete_ = false;
    complete_points_abs();


}

Optimization::Case* PolyModel::CaseFromPoint(Eigen::VectorXd point, Optimization::Case *prototype) {

    //std::cout<<" Create new case from point "<< point<<std::endl;
    // In order for case to exist outside poly_model, we use the new operator
    Optimization::Case *new_case = new Optimization::Case(prototype);
    new_case->SetRealVarValues(point);
    new_case->set_objective_function_value(std::numeric_limits<double>::max());
    //objective value must be changed to not evalated, i.e. MAXLIMITSDOUBLE
    return new_case;
}

Eigen::VectorXd PolyModel::find_new_point(Polynomial basis_function) {

    int dimension = basis_function.return_dimension();
    Eigen::VectorXd coeffs = basis_function.return_coeffs();
    Eigen::VectorXd x0, x1, x2, x3, x4;
    x0 = x1 = x2 = x3 = x4 = Eigen::VectorXd::Zero(dimension);

    // Find largest monomial coefficient (excluding constant term which has already been assigned to first point)
    double max = 0.0;
    int max_coeff = -1;
    for (int i = 1; i < coeffs.size(); ++i) {
        if(fabs(coeffs(i)) > max) {
            max = fabs(coeffs(i));
            max_coeff = i;
        }
    }

    if(max_coeff==-1){
        std::cout << "Good point alg, Problem: all coefficients are zero, should never happen" << std::endl;
        // LETS JUST PRINT A FEW OF THE COEFFS
        for(int i=0; i<3; i++){
            std::cout << coeffs[i] << std::endl;
        }
    }
    if(max_coeff<=dimension){
        // The largest coefficient is from a linear term
        x1(max_coeff-1) = 1; //subract 1 maybe if change polynomial form
        x2 = -x1;
    }
    else if(max_coeff<=2*dimension){
        // Largest coefficient is quadratic monomial
        x1(max_coeff-dimension-1) = 1;
        x2 = -x1;
    }
    else{
        // Mixed term quadratic is largest
        // There is probably a smarter way to do this... oh whatever I'm lazy
        int l,m = -1;
        int coeff_dummy = 2*dimension+1;

        for(int i=0; i<dimension-1; i++){
            for (int j=i+1; j<dimension; j++) {
                if (max_coeff == coeff_dummy) {
                    l = i;
                    m = j;
                    goto endloop;
                }
                coeff_dummy = coeff_dummy+1;
            }
        }
        endloop:

        x1(l) =  1.0/sqrt(2);
        x1(m) = -1.0/sqrt(2);
        x2 = -x1;
        for(int i=0; i<dimension; i++){
            x3(i) = fabs(x1(i));
            x4(i) = -fabs(x1(i));
        }
    }
    Eigen::VectorXd best_point;
    double best_value = 0.0;
    QList<Eigen::VectorXd> points;
    points.append(x0);
    points.append(x1);
    points.append(x2);
    points.append(x3);
    points.append(x4);

    // Determine which of the 5 points is the best one
    for(int i=0; i<5; i++){
        if(fabs(basis_function.evaluate(points.at(i)))>=best_value){
            best_point = points.at(i);
            best_value = fabs(basis_function.evaluate(points.at(i)));
        }
    }

    return best_point;
}
//complete interpolation points scaled to a ball of radius 1 around origin

void PolyModel::complete_points_abs() {
    Eigen::VectorXd centre_point = points_.at(0);
    // points_abs=scaled interpolation point.
    // points_abs.at(0) is always zero
    points_abs.append(points_.at(0) - centre_point);
    //std::cout<<" the first element in points_abs is "<<points_abs.at(0)<<std::endl;

    int n_Polynomials = basis_.length();
    int n_points = points_.length();
    QList<Polynomial> temp_basis = get_basis();

    for (int i = 0; i < n_Polynomials; i++) {
        Polynomial cur_pol = temp_basis.at(i);
        if(i==0){
            /*std::cout << "we don't need to find new point, basis polynomial, i = " << i << std::endl;
            std::cout << " This point is "<<points_abs.at(i)<<std::endl;*/
        }
        else{
        /*std::cout << "we need to find new point, basis polynomial, i = " << i << std::endl;*/
        Polynomial temp_poly_here = temp_basis.at(i);

        // Append new point and swap it to current position
        points_abs.append(find_new_point(temp_poly_here));
        //points_abs.swap(i, points_abs.length() - 1);
        /*std::cout<<"this point is "<< points_abs.at(i)<<std::endl;*/

        Polynomial temp_i = temp_basis.at(i);
        auto temp_point = points_abs.at(i);

        for (int j = i + 1; j < n_Polynomials; j++) {
            Polynomial uj = temp_basis.at(j);
            Polynomial ui = temp_basis.at(i);
            double ratio;
            //TODO: TESTING DIVISION BY ZERO
            if(ui.evaluate(temp_point)==0 && j==i+1){
                std::cout << "Division by zero because U_i(y_i) = 0" << std::endl;}
            // END TESTING
            if(uj.evaluate(temp_point) == ui.evaluate(temp_point)){
                ratio = 1.0;
            }
            else{ratio = uj.evaluate(temp_point) / ui.evaluate(temp_point);}
            ui.multiply(-1.0 * ratio);
            uj.add(ui);
            temp_basis[j] = uj;
        }}

    }

}

// scale points back, complete interpolation points
void PolyModel::complete_points() {
    int n_Polynomials = basis_.length();
    //std::cout<<"n_Plynomials is "<<n_Polynomials<<std::endl;
    Eigen::VectorXd centre_point = points_.at(0);
    /*std::cout<< "for i= 0" <<std::endl;
    std::cout<< "This interpolation point is:\n  " <<points_.at(0)<<std::endl;*/
    // Scale points back
    for (int i = 1; i < n_Polynomials; ++i) {
        points_.append(centre_point + radius_*points_abs.at(i));
      /*  std::cout<< "for i= " <<i<<std::endl;
        std::cout<< "This interpolation point is:\n " <<points_.at(points_.size()-1)<<std::endl;*/

        //creat case from new interpolation point

        cases_.append(CaseFromPoint(centre_point + radius_*points_abs.at(i), cases_.at(0)));
        //std::cout<< "  for i=  " << i << " point of cases_.at(i) is "<< cases_.at(i)->GetRealVarVector()<<std::endl;

        // Append case to list of unevaluated cases
         cases_not_eval_.append(cases_.at(i));
        //std::cout<< " send point to cases_not_eval "<< cases_not_eval_.at(cases_not_eval_.size()-1)->GetRealVarVector()<<std::endl;

        //Eigen::VectorXd point=cases_not_eval_.at(end)->GetRealVarVector();
         }
    needs_evals_ = true;
    needs_set_of_points_ = false;
    std::cout<< "Interpolation points are completed:\n  " <<std::endl;
}

void PolyModel::calculate_model_coeffs() {
    std:: cout <<"Calculate model coefficent " << std::endl;
    if(!needs_evals_ && !needs_set_of_points_){
        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(basis_.length() ,basis_.length());
        Eigen::VectorXd y = Eigen::VectorXd::Zero(basis_.length());


        // Build Matrix M from basis and functions evaluations
        for (int i = 0; i < basis_.length(); ++i) {
            for (int j = 0; j < basis_.length(); ++j) {
                Polynomial current_basis = basis_.at(j);
                M(i,j) = current_basis.evaluate(points_.at(i));
            }
            y(i) = cases_.at(i)->objective_function_value();
        }

        Eigen::VectorXd alpha = M.inverse()*y;
        model_coeffs_ = alpha;
        is_model_complete_= true;
        std::cout << "Current Model coefficient is\n" <<model_coeffs_<<std::endl;

    }
    else{std::cout << "Model_coefficient alg: Either needs evaluations or set of points not finished yet" << std::endl;}
}

// Find optimization step along Gradient
Eigen::VectorXd PolyModel::optimizationStep_CP() {
    std::cout << "Calculate OptimizationStep" <<std::endl;
   Eigen::VectorXd coeffs= get_model_coeffs();
    Polynomial Poly = Polynomial(dimension_, coeffs);
    grad=Poly.evaluateGradient(center_);
    //Find Hessian matrix first
    Poly.Hessian();
    //std::cout<<"Find Cauchy point:"<<std::endl;
    optimization_step_CP= Poly.Cauchy_Point(center_,radius_,grad); // Optimizationstep at subregion boundary based on G
    //std::cout << " OptimizationStep(Gradient Descent) is" <<optimization_step_CP <<std::endl;
    return optimization_step_CP;
}
// Return the intersection of the segment
// connecting the Cauchy point to the Newton point with the
// trust region boundary
Eigen::VectorXd PolyModel::optimizationStep_SDL() {
    std::cout << "Calculate OptimizationStep" <<std::endl;
    Eigen::VectorXd coeffs= get_model_coeffs();
    Polynomial Poly = Polynomial(dimension_, coeffs);
    grad=Poly.evaluateGradient(center_);
    //Find Hessian matrix first
    Poly.Hessian();
   optimization_step_SDL= Poly.Dogleg_step(center_,radius_,grad); // Optimizationstep at subregion boundary based on G
    //std::cout << " OptimizationStep(Gradient Descent) is" <<optimization_step_CP <<std::endl;

    return optimization_step_SDL;

}
Eigen::VectorXd PolyModel::Gradient(){
    Eigen::VectorXd coeffs= get_model_coeffs();
    Polynomial Poly = Polynomial(dimension_, coeffs);
    grad=Poly.evaluateGradient(center_);
    return grad;
}



void PolyModel::addBaseCase( Optimization::Case* BaseCase) {
   // Optimization::Case *newBaseCase=CaseFromPoint(NewCenterPoint,cases_.at(0));//get new BaseCase from NewCenterPoint
    //clear points_ and cases_
    cases_ = QList<Optimization::Case*>();
    points_= QList<Eigen::VectorXd> ();

    // add newBaseCase to cases_and cases_not_eval
    cases_.append(BaseCase);
    //cases_not_eval_.append(newBaseCase);
    points_.append(BaseCase->GetRealVarVector());
    center_ = points_.at(0);
    std::cout << "The new center point is " <<center_<<std::endl;

    is_model_complete_ = false;
    needs_set_of_points_ = true;
    needs_evals_= true;

}


double PolyModel::obejctive_function_value_model(Eigen::VectorXd point) {
    Polynomial polynomial=Polynomial(dimension_,get_model_coeffs());
    double objective_function_value_model=polynomial.evaluate(point);
    return objective_function_value_model;

}