//
// Created by cutie on 07.10.16.
//

#include "polynomial.h"

Polynomial::Polynomial(int dimension, Eigen::VectorXd coeffs) {
    dimension_ = dimension;
    coeffs_ = coeffs;
    no_elemts_ = (dimension+1)*(dimension+2)/2;
}

double Polynomial::evaluate(Eigen::VectorXd point) {
    /* Assume quadratic monomial basis as presented
 * in Caio Giuliani's summary, but with a slight
 * change in ordering. Here the constant term is
 * the first one, then the linear terms, then
 * the quadratic terms of one variable, and finally
 * the mixed terms. i.e. (1, x_1, x_2,..., 0.5*x_N,
 * 0.5*x_1^2, x_2^2,..., 0.5*x_N^2, x_1*x_2, x_1*x_3,
 * ..., x_(N-1)*x_N )
 */

    // CHECK POLYNOMIAL AND POINT DIMENSIONS SAME
    double sum = coeffs_(0);
    for(int i=0; i<dimension_; i++){
        // Linear term and quadratic single term
        sum = sum + coeffs_(i+1)*point(i);
        sum = sum + 0.5*coeffs_(dimension_+i+1)*point(i)*point(i);
    }

    // iter is just an iterator to keep track of where we are in the
    // coefficient vector. Simpler than writing as a function of i and j
    int iter = 2*dimension_+1;

    for(int i=0; i<dimension_-1; i++){
        for(int j=i+1; j<dimension_; j++){
            sum = sum + coeffs_(iter)*point(i)*point(j);
            iter = iter +1;
        }
    }

    return sum;
}

Eigen::VectorXd Polynomial::evaluateGradient(Eigen::VectorXd point) {

    /* We use the taylor expansion of a function and the
     * fact that our polynomial's third (and higher) derivative
     * is zero, so the truncation error of our approximaton is
     * zero as well and thus the approximation is actually an
     * exact evaluation.
     */

    Eigen::VectorXd grad(dimension_);
    for (int i = 0; i < dimension_; ++i) {

        // We could choose any length but unit length is simple
        Eigen::VectorXd unit = Eigen::VectorXd::Zero(dimension_);
        unit(i) = 1;

        // Central difference approximation
        double grad_value_i = 0.5*(evaluate(point+unit)-evaluate(point-unit));
        grad(i) = grad_value_i;
    }
    return grad;

}


Eigen::MatrixXd  Polynomial::Hessian(){
    std::cout<<"Find Hessian Matrix"<<std::endl;
    /* Hessian matrix is a square matrix of second-order
     * partial derivative of a scalar-valued function.
     * for quadratic poly model, Hessian matrix only depend on cofficients
 */
    Eigen::MatrixXd B=Eigen::MatrixXd::Zero(dimension_,dimension_);

    for(int i=0; i<dimension_; ++i){
            B(i,i)=coeffs_(dimension_+i+1);
    }

    for (int i=0;i<dimension_; ++i){
        int k=0;
        for (int j=i+1;j<dimension_; ++j) {
                int a = ((dimension_ + 1) + (dimension_ - i)) * (i + 2) / 2 + k;;
                B(i,j) = coeffs_(a);
                B(j,i)=B(i,j);
                 k=k+1;

        }
        }

    std::cout<< "Here is the Hessian matrix B :\n"<< B <<std::endl;
    Hessian_Matrix=B;

 return Hessian_Matrix;
}

Eigen::VectorXd Polynomial::Cauchy_Point(Eigen::VectorXd points)
{
   // Eigen::MatrixXd Hessian_Matrix=Hessian();
    Eigen::MatrixXd Inverse_Matrix=Hessian_Matrix.inverse();
    std::cout<<"Find inverse Hessian Matrix:\n"<<Inverse_Matrix<<std::endl;
    Eigen::VectorXd grad=evaluateGradient(points);
    Eigen::VectorXd Newton_Point=-Inverse_Matrix*grad;
    std::cout<<"The Quasi-Newton Point is:\n"<<Newton_Point<<std::endl;
    std::cout<<"length of Newton Point is:\n"<<Newton_Point.norm()<<std::endl;
    std::cout<<"objective function value of current Point is:\n"<<evaluate(points)<<std::endl;
    std::cout<<"objective function value of Newton Point is:\n"<<evaluate(Newton_Point)<<std::endl;
    return Newton_Point;

}

Eigen::VectorXd Polynomial:: Newton_Point(Eigen::VectorXd points) {

    //Eigen::MatrixXd Hessian_Matrix=Hessian();
    Eigen::MatrixXd Inverse_Matrix=Hessian_Matrix.inverse();
    std::cout<<"Find inverse Hessian Matrix:\n"<<Inverse_Matrix<<std::endl;
    Eigen::VectorXd grad=evaluateGradient(points);
    Eigen::VectorXd Newton_Point=-Inverse_Matrix*grad;
    std::cout<<"The Quasi-Newton Point is:\n"<<Newton_Point<<std::endl;
    std::cout<<"length of Newton Point is:\n"<<Newton_Point.norm()<<std::endl;
    std::cout<<"objective function value of current Point is:\n"<<evaluate(points)<<std::endl;
    std::cout<<"objective function value of Newton Point is:\n"<<evaluate(Newton_Point)<<std::endl;
    return Newton_Point;

}


void Polynomial::add(Polynomial poly) {

    if(poly.return_no_elements()!= no_elemts_){
        std::cout << "polynomials have different dimensions" << std::endl; //TODO Throw exception
    }

    else{
        coeffs_ = coeffs_ + poly.coeffs_;
    }
}

void Polynomial::multiply(double k) {
    coeffs_ = k*coeffs_;
}
