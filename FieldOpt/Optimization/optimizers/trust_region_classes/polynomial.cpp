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

Eigen::VectorXd Polynomial::Cauchy_Point(Eigen::VectorXd points, double radius_)
{      Eigen::VectorXd grad=evaluateGradient(points);

    double gBg=grad.transpose()*Hessian_Matrix*grad;
    //std::cout<<"gBg is: "<<gBg<<std::endl;
    double tau;
   if(gBg <=0) {
       std::cout<<"gbg <0 "<<std::endl;
       tau=1;
   }
   else {
       double length=grad.norm()*grad.norm()*grad.norm()/(radius_*gBg);
       if (length>=1){
           std::cout<<"The Cauchy Point lies outside the trust region, then take a step to the boundary along steepest descent direction:"<<std::endl;
           tau=1;
       }
       else{
           std::cout<<"The Cauchy Point lies inside the trust region,then take Cauchy Point as optimization step size:"<<std::endl;
           tau=length;
       }

   }
    Eigen::VectorXd Cauchy_Point=-tau*radius_*grad/grad.norm();
    std::cout<<"tau is: "<<tau<<std::endl;
    std::cout<<"The Cauchy_Point is:\n"<<Cauchy_Point<<std::endl;
    std::cout<<"length of Cauchy Point is:\n"<<Cauchy_Point.norm()<<std::endl;
    //std::cout<<"objective function value of current Point is:\n"<<evaluate(points)<<std::endl;
    //std::cout<<"objective function value of Cauchy Point is:\n"<<evaluate(Cauchy_Point)<<std::endl;
    return Cauchy_Point;

}


Eigen::VectorXd Polynomial:: Newton_Point(Eigen::VectorXd points, double radius_) {

    //Eigen::MatrixXd Hessian_Matrix=Hessian();
    Eigen::MatrixXd Inverse_Matrix=Hessian_Matrix.inverse();
    //std::cout<<"Find inverse Hessian Matrix:\n"<<Inverse_Matrix<<std::endl;
    Eigen::VectorXd grad=evaluateGradient(points);
    Eigen::VectorXd Newton_Point=-Inverse_Matrix*grad;
    std::cout<<"The Quasi-Newton Point is:\n"<<Newton_Point<<std::endl;
    std::cout<<"length of Newton Point is:\n"<<Newton_Point.norm()<<std::endl;
    //std::cout<<"objective function value of current Point is:\n"<<evaluate(points)<<std::endl;
    //std::cout<<"objective function value of Newton Point is:\n"<<evaluate(Newton_Point)<<std::endl;
    if (Newton_Point.norm()>radius_)
    {
        std::cout<<"The Newton Point lies outside the trust region"<<std::endl;
    }
    else{
        std::cout<<"The Newton Point lies inside the trust region"<<std::endl;
    }
    return Newton_Point;

}

// Newton point lies inside trustregion, cauchy point lies outside trustregion

Eigen::VectorXd  Polynomial::Dogleg_step(Eigen::VectorXd points,double radius_)

{
    Eigen::VectorXd Dogleg_step;
    std::cout<<"Find Cauchy point:\n"<<std::endl;
    Eigen::VectorXd d_cp= Cauchy_Point(points,radius_);
    std::cout<<"Find Quastion-Newton point:\n"<<std::endl;
    Eigen::VectorXd d_nw=Newton_Point(points, radius_);
    if(d_cp.norm()>=radius_)
    {
        Dogleg_step=d_cp;
        std::cout << " Cauchy point lies on the boundary, then take Cauchy point as optimization step:\n" << Dogleg_step <<std::endl;
    }
    else if (d_nw.norm()<=radius_){
        Dogleg_step=d_nw;
        std::cout << " Both Cauchy point and Quastion-Newton point lie inside trust region, then take Newton point as optimization step:\n" << Dogleg_step <<std::endl;

    }
    else
    {
        std::cout << "Cauchy point lies inside trust region, Quastion-Newton point lies outside trust region, then take Dogleg optimization step:\n" << Dogleg_step <<std::endl;
         std::cout<<"Find tau such that || d_cp + tau*(d_nw - d_cp) || = radius :\n"<<std::endl;
        // solve a single quadratic equation
        double tau;
        Eigen::VectorXd diff=d_nw-d_cp;
        double a=diff.dot(diff);
        double b=2*d_cp.dot(diff);
        double c=d_cp.dot(d_cp)-radius_*radius_;
        tau=(-b+sqrt(b*b-4*a*c))/(2*a);
         std::cout<<"tau is:\n"<<tau<<std::endl;
        Dogleg_step=d_cp+tau*(d_nw-d_cp);
        std::cout << "Dogleg optimization step is:\n" << Dogleg_step <<std::endl;
        std::cout<<"length of Dogleg point is:\n"<<Dogleg_step.norm()<<std::endl;
    }
    return Dogleg_step;

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
