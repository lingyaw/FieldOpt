//
// Created by cutie on 07.10.16.
//

#ifndef FIELDOPT_POLY_MODEL_H
#define FIELDOPT_POLY_MODEL_H

#include "Optimization/case.h"
#include "polynomial.h"
#include <QList>
#include <iostream>

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
    QList<Optimization::Case*> cases_; //!< All cases, might or might not be evaluated
    QList<Optimization::Case*> cases_not_eval_; //!< Cases that need to be sent to case_handler

    // Bools to help trust_region_search class determine next step
    bool needs_evals_; //
    bool needs_set_of_points_;
    bool is_model_complete_; //!< Bool that is set true whenever a model has been build for current center point and radius

    Eigen::VectorXd center_;
    double radius_;
    int dimension_;
    QList<Polynomial> basis_; //!< Monomial basis of model, usually quadratic
    Eigen::VectorXd model_coeffs_; //!< The coefficients of the model using basis
    Eigen::VectorXd optimization_step_NG;
    Eigen::VectorXd optimization_step_SDL;

    /*!
     * @brief As described by A. Conn, finds a 'good point' for the
     * scaled trust region. This is a copy of C. Giuliani's Matlab code
     * and will always find a point such that abs(basis_function(x)) > 0.24
     * @param Double k
     * @return A good point
     */
    Eigen::VectorXd find_new_point(Polynomial basis_function);

public:
    /*!
     * @brief PolyModel constructor.
     */
    PolyModel(Optimization::Case* base_case, double radius);

    /*!
     * @brief Empty default PolyModel constructor so other classes can contain a PolyModel as private variable
     */
    PolyModel() {};

    /*!
     * \brief create a Case from a list of variables and a Case prototype
     *
     * Creates a Case type object from a Case prototype (i.e. a case with the same
     * number of variables but where the variable values have been altered.
     * \return A Case generated from a Eigen::VectorXd point
     */
    static Optimization::Case* CaseFromPoint(Eigen::VectorXd point, Optimization::Case *prototype);
    //Optimization::Case*   find_NewBaseCase();

    QList<Eigen::VectorXd> get_points() {
        return points_;
    };

    QList<Optimization::Case*> get_cases_not_eval() {
        return cases_not_eval_;
    };

    void ClearCasesNotEval() {
        cases_not_eval_ = QList<Optimization::Case*>();
        needs_evals_ = false;
    }

    bool get_needs_evals() {
        return needs_evals_;
    }

    double get_radius() {
        return radius_;
    };

    QList<Polynomial> get_basis() {
        return basis_;
    };

    Eigen::VectorXd get_model_coeffs() {
       // calculate_model_coeffs();
        return model_coeffs_;
    };

    //Eigen::VectorXd get_optimizationstep() {
       // optimizationStep();
      //  return optimization_step;
   // };


    bool ModelNeedsEvals() {
        return needs_evals_;
    }

    bool ModelNeedsSetOfPoints() {
        return needs_set_of_points_;
    }

    bool isModelReady() {
        return is_model_complete_;
    }

    void setRadius(double k){
        radius_ = k;
    };

    /*!
     * @brief Complete set of interpolation points
     * using Algorithm 5 as described in paper
     * by C. Giuliani
     */
    void complete_points();

    /*!
     * @brief Calculate coefficients of quadratic model
     * of the trust region using a complete and
     * well poised set of points
     */
    void calculate_model_coeffs();

    void set_evaluations_complete(){
        needs_evals_ = false;
    }
    void set_model_complete(){
        is_model_complete_ = true;
    }
    void addCenterPoint(Eigen::VectorXd newCenterPoint);


    Eigen::VectorXd optimizationStep_NG();
    Eigen::VectorXd optimizationStep_SDL();

    Eigen::VectorXd get_centerpoint(){
        return center_;
    }

};

#endif //FIELDOPT_POLY_MODEL_H
