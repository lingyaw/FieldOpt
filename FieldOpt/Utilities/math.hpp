/// This file contains helper methods for working with maths.
#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include <QList>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <boost/random.hpp>

/*!
 * @brief Calculate the average value of the items in the list. The returned value will always be a double.
 * @param list
 * @return The average value of the elements in the list as a double.
 */
template<typename T>
inline double calc_average(const QList<T> list) {
    assert(!list.empty());
    return std::accumulate(list.begin(), list.end(), 0.0) / list.size();
}

/*!
 * @brief Calculate the median of a list. The returned value will have the same type as the elements in the list.
 * If the values are integers and the list has an even number of elements, the result will be floored.
 * @param list
 * @return
 */
template<typename T>
inline T calc_median(QList<T> list) {
    assert(!list.empty());
    std::sort(list.begin(), list.end());
    size_t size = list.size();
    if (size % 2 == 0)
        return (list[size/2 - 1] + list[size/2])/2;
    else
        return list[size/2];
}

/*!
 * @brief The Polynomial class is an implementation to describe second
 * order polynomials with the natural monomial basis of second order polynomials.
 */
class Polynomial {
private:
    int dimension_;
    int no_elemts_;
    Eigen::VectorXd coeffs_;

public:
    /*!
     * @brief Polynomial constructor. Sets coefficients, dimension and calculates number of elements.
     * @param dimension of polynomial.
     * @param An Eigen::VectorXd containing coefficients of polynomial.
     * @return
     */
    Polynomial(int dimension, Eigen::VectorXd coeffs) {
        dimension_ = dimension;
        coeffs_ = coeffs;
        no_elemts_ = (dimension+1)*(dimension+2)/2;
    };

    /*!
     * @brief Empty default Polynomial constructor so other classes can contain a Polynomial as private variable
     */
    Polynomial();

    Eigen::VectorXd return_coeffs() {
        return coeffs_;
    };
    int return_dimension() {
        return dimension_;
    };
    int return_no_elements() {
        return no_elemts_;
    };

    /*!
     * @brief Evaluates polynomial in a given input point.
     * @param evaluation point.
     * @return Polynomial value.
     */
    double evaluate(Eigen::VectorXd point) {
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
    };

    /*!
     * @brief Evaluates gradient of model in a given input point.
     * @param gradient evaluation point.
     * @return gradient value.
     */
    Eigen::VectorXd evaluateGradient(Eigen::VectorXd point) {

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

    };

    /*!
     * @brief Adds (the coefficients of) a polynomial to *this polynomial
     * @param Polynomial
     */
    void add(Polynomial poly) {

        if(poly.return_no_elements()!= no_elemts_){
            std::cout << "polynomials have different dimensions" << std::endl;
        }

        else{
            coeffs_ = coeffs_ + poly.coeffs_;
        }

    };

    /*!
     * @brief Multiplies all coefficients of current polynomial with an input double
     * @param Double k
     */
    void multiply(double k) {
        coeffs_ = k*coeffs_;
    };

};

/*!
 * @brief The PolyModel class is a model for describing a
 * function with polynomials within a given radius of a
 * center point (i.e. trust region). Before the model is
 * computed we must provide a function that is to be
 * approximated and as many function evaluations as we
 * wish. The PolyModel class can then calculate the needed
 * points in order to make the polynomial fitting unique.
 */
class PolyModel {

private:
    // Member variables
    QList<Eigen::VectorXd> points_;
    QList<int> points_evaluated_;
    bool needs_evals_;
    QList<double> fvalues_;
    Eigen::VectorXd center_;
    double radius_;
    int dimension_; // Monomial basis of model, usually quadratic
    QList<Polynomial> basis_; // The coefficients of the model using basis
    Eigen::VectorXd model_coeffs_; // Private methods/sub-methods
    bool is_model_complete_; // Bool that is set true whenever a model has been build for current center point and radius

   /*!
    * @brief As described by A. Conn, finds a 'good point' for the
    * scaled trust region. This is a copy of C. Giuliani's Matlab code
    * @param Double k
    * @return A good point
    */
    Eigen::VectorXd find_new_point(Polynomial poly) {

    int dimension = poly.return_dimension();
    Eigen::VectorXd coeffs = poly.return_coeffs();
    Eigen::VectorXd x0, x1, x2, x3, x4;
    x0 = x1 = x2 = x3 = x4 = Eigen::VectorXd::Zero(dimension);

    // Find largest monomial coefficient (excluding constant term which has already been assigned to first point)
    double max = 0.0;
    int max_coeff = -1;
    for (int i = 1; i < poly.return_no_elements(); ++i) {
        if(fabs(coeffs(i)) > max) {
            max = fabs(coeffs(i));
            max_coeff = i;
        }
    }
    if(max_coeff==-1){ std::cout << "Good point alg, Problem: all coefficients are zero, should never happen" << std::endl;}
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
            for (int j=1; j<dimension; j++) {
                if (max_coeff == coeff_dummy) {
                    l = i;
                    m = j;
                    break;
                }
                coeff_dummy = coeff_dummy+1;
            }
        }

        x1(l) =  1.0/sqrt(2);
        x1(m) = -1.0/sqrt(2);
        x2 = -x1;
        for(int i=0; i<dimension; i++){
            x3(i) = fabs(x1(i));
            x4(i) = -fabs(x1(i));
        }
    }
    Eigen::VectorXd best_point;
    double best_value = 0;
    QList<Eigen::VectorXd> points;
    points.append(x0);
    points.append(x1);
    points.append(x2);
    points.append(x3);
    points.append(x4);

    // Determine which of the 5 points is the best one
    for(int i=0; i<5; i++){
        if(fabs(poly.evaluate(points.at(i)))>=best_value){
            best_point = points.at(i);
            best_value = fabs(poly.evaluate(points.at(i)));
        }
    }

    return best_point;
    };

public:
    /*!
     * @brief PolyModel constructor.
     */
    PolyModel(QList<Eigen::VectorXd> points, QList<double> fvalues, double radius, int dimension) {
        points_ = points;
        fvalues_ = fvalues;
        // Add 1s to signal that points have been evaluated
        for(int i=0; i<fvalues.length(); i++){points_evaluated_.append(1);}
        needs_evals_ = false; // All points have been evaluated
        center_ = points.at(0);
        radius_ = radius;
        dimension_ = dimension;
        is_model_complete_ = false;
        QList<Polynomial> basis;
        for (int i = 0; i < (dimension+1)*(dimension+2)/2; ++i) {
            Eigen::VectorXd temp_vec = Eigen::VectorXd::Zero((dimension+1)*(dimension+2)/2);
            temp_vec(i) = 1;
            Polynomial temp_poly = Polynomial(dimension, temp_vec);
            basis.append(temp_poly);
        }
        basis_ = basis;
    };

    /*!
     * @brief Empty default PolyModel constructor so other classes can contain a PolyModel as private variable
     */
    PolyModel();

    QList<Eigen::VectorXd> get_points() {
        return points_;
    };

    QList<int> get_points_evaluated() {
        return points_evaluated_;
    };

    bool get_needs_evals() {
        return needs_evals_;
    }

    QList<double> get_fvalues() {
        return fvalues_;
    };

    double get_radius() {
        return radius_;
    };

    QList<Polynomial> get_basis() {
        return basis_;
    };

    Eigen::VectorXd get_model_coeffs() {
        return model_coeffs_;
    };

    bool isModelComplete() {
        return is_model_complete_;
    }

    /*!
     * @brief Complete set of interpolation points
     * using Algorithm 5 as described in paper
     * by C. Giuliani
     */
    void complete_points(){
        // Complete set of interpolation points so
        // that the set is well-poised

        /* Pivot element tolerance. Note that we can always find
         * an element inside the ball of radius one in order to
         * get as close as we want to 0.25. It might however be clever
         * to lower this tolerance in order to preserve as many stored
         * function evaluations as possible.
         */
        double tol_pivot = 0.24;
        Eigen::VectorXd centre_point = points_.at(0);

        /* scaling points to a ball of radius 1
         * and center at center_ (first point of
         * the points list)
         */

        QList<Eigen::VectorXd> points_abs;
        for (int i = 0; i < points_.length(); ++i) {
            points_abs.append((1/radius_)*(get_points().at(i) - centre_point));
        }

        int n_Polynomials = basis_.length();
        int n_points = points_.length();
        QList<Polynomial> temp_basis = get_basis();

        for(int i=0; i<n_Polynomials; i++){
            Polynomial cur_pol = temp_basis.at(i);
            double max_abs = 0.0;
            int max_abs_ind = -1;
            for (int j = i; j < n_points; ++j) {
                /* If new max value, and within a distance radius_ of center point,
                 * accept this point as the currently best point.
                 * Note that we have scaled the points so we only need to
                 * check that the norm of the current points is <=1
                 */
                if(fabs(cur_pol.evaluate(points_abs.at(j)))>max_abs && points_abs.at(i).norm()<=1){
                    max_abs = fabs(cur_pol.evaluate(points_abs.at(j)));
                    max_abs_ind = j;
                }
            }

            /* If evaluation in pivot element is greater than threshold,
             * switch elements
             * and its associated function evaluations.
             */
            if(max_abs>tol_pivot) {
                //YES sufficient pivot element aka. good point
                points_abs.swap(i,max_abs_ind);
                fvalues_.swap(i,max_abs_ind);
                points_evaluated_.swap(i,max_abs_ind);
            }
            else{
                //NO sufficient pivot element aka. good point
                //Find new point using alg proposed by Conn
                Polynomial temp_poly_here = temp_basis.at(i);
                // Append new point and swap it to current position
                points_abs.append(find_new_point(temp_poly_here));
                points_abs.swap(i,points_abs.length()-1);
                // Append 0 to signal that point is not yet evaluated
                points_evaluated_.append(0);
                points_evaluated_.swap(i,points_abs.length()-1);
                // Append dummy value and change needs_evals_ to true;
                fvalues_.append(silly_function(points_abs.at(i)));
                fvalues_.swap(i,points_abs.length()-1);
                needs_evals_ = true;


            }

            Polynomial temp_i = temp_basis.at(i);
            auto temp_point = points_abs.at(i);


            for (int j = i+1; j <n_points ; j++) {
                Polynomial uj= temp_basis.at(j);
                Polynomial ui = temp_basis.at(i);
                double temp_ratio = uj.evaluate(temp_point)/ui.evaluate(temp_point);
                ui.multiply(-1.0*temp_ratio);
                uj.add(ui);
                temp_basis[j] = uj;
            }
        }

        // Scale points back
        QList<Eigen::VectorXd> points_scaled;
        for (int i = 0; i < n_Polynomials; ++i) {
            points_scaled.append(centre_point + radius_*points_abs.at(i));
        }

        points_ = points_scaled;
    };

    /*!
     * @brief Calculate coefficients of quadratic model
     * of the trust region using a complete and
     * well poised set of points
     */
    void calculate_model_coeffs() {
        if(needs_evals_){std::cout << "there are unfinished evaluations, model will be wrong" << std::endl;}

        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(basis_.length() ,basis_.length());
        Eigen::VectorXd y = Eigen::VectorXd::Zero(basis_.length());


        // Build Matrix M from basis and functions evaluations
        for (int i = 0; i < basis_.length(); ++i) {
            for (int j = 0; j < basis_.length(); ++j) {
                Polynomial current_basis = basis_.at(j);
                M(i,j) = current_basis.evaluate(points_.at(i));
            }
            //y(i) = fvalues_.at(i);
            y(i) = silly_function(points_.at(i));
        }

        Eigen::VectorXd alpha = M.inverse()*y;
        model_coeffs_ = alpha;
        is_model_complete_ = true;

    };

    /*!
    * @brief Returns a list of indeces of all non-evaluated cases
    * @return List of indeces of non-evaluated cases
    */
    QList<int> missing_evaluations() {
        QList<int> missing;
        for(int i=0; i<points_evaluated_.length(); i++){
            if(points_evaluated_.at(i)==0){missing.append(i);}
        }
        return missing;
    }

    void add_point(Eigen::VectorXd point){
        points_.append(point);
        points_evaluated_.append(0);
        needs_evals_ = true;
    }

    void set_fvalue(double fvalue, int i){
        fvalues_[i]=fvalue;
        points_evaluated_[i]=1;
    }

    void set_evaluations_complete(){
        needs_evals_ = false;
    }


    // Silly function for testing:
    double silly_function(Eigen::VectorXd x) {
        return 3+ 4*x(0) + 3*x(1) + x(0)*x(0) + 5*x(1)*x(1) -1*x(0)*x(1);
    };

};

 * @brief Linspace function. Create a list containing all elements in the range from start to (but not including) end
 * with step between each element.
 *
 * Examples: range(1, 10, 1) = [1, 2, 3, 4, 5, 6, 7, 8, 9]
 *           range(0, 10, 3) = [0, 3, 6, 9]
 * @tparam T A number type, e.g. int or double.
 * @param start First number in range.
 * @param end End of range.
 * @param step Step length.
 * @return
 */
template <typename T>
inline std::vector<T> range(T start, T end, T step) {
    int length = abs((end - start) / step);
    auto ran = std::vector<T>(length);
    ran[0] = start;
    for (int i = 1; i < length; ++i) {
        ran[i] = ran[i-1] + step;
    }
    return ran;
}

/*!
 * @brief Get a random integer in the range [lower .. upper]
 * @param lower The lowest possible int.
 * @param upper The highest possible int.
 * @return A random integer.
 */
inline int random_integer(const int lower, const int upper) {
    boost::random::mt19937 gen;
    boost::random::uniform_int_distribution<> dist(lower, upper);
    return dist(gen);
}

#endif // MATH_FUNCTIONS_H
