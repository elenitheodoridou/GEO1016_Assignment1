/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{
    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
                 "\t    - calibration_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
                 "\t    - compute the SVD of a matrix;\n"
                 "\t    - compute the inverse of a matrix;\n"
                 "\t    - compute the transpose of a matrix.\n\n"
                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
                 "\tfunctions are provided in the following files:\n"
                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    std::cout << "\n[Liangliang]:\n"
                 "\tThe input parameters of this function are:\n"
                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
                 "\tparameters are returned by the following variables:\n"
                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t:         a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // Input validation
    if (points_2d.size() == points_3d.size() && points_2d.size() >= 6) {
        std::cout << "Input is valid" << std::endl;

        // Computing the P matrix (so P * m = 0)
        int number_of_equations = points_3d.size() * 2;
        Matrix matrix_P(number_of_equations, 12, 0.0);

        int k = 0;
        for (int i = 0; i < points_3d.size(); i++) {
            std::vector<double> vector_ax = {points_3d[i][0], points_3d[i][1], points_3d[i][2], 1, 0, 0, 0, 0, -points_2d[i][0]*points_3d[i][0], -points_2d[i][0]*points_3d[i][1], -points_2d[i][0]*points_3d[i][2], -points_2d[i][0]};
            std::vector<double> vector_ay = {0, 0, 0, 0, points_3d[i][0], points_3d[i][1], points_3d[i][2], 1, -points_2d[i][1]*points_3d[i][0], -points_2d[i][1]*points_3d[i][1], -points_2d[i][1]*points_3d[i][2], -points_2d[i][1]};
            if (i == 0) {
                matrix_P.set_row(i, vector_ax);
            } else {
                matrix_P.set_row(i+k, vector_ax);
            }
            k++;
            matrix_P.set_row(i+k, vector_ay);
        }

        std::cout << "Matrix P \n" << matrix_P << std::endl;

        // Solve matrix M (the whole projection matrix, i.e. M = K*[R, t]) using SVD decomposition
        // SVD decomposition
        Matrix matrix_U(number_of_equations, number_of_equations, 0.0);
        Matrix matrix_S(number_of_equations, 12, 0.0);
        Matrix matrix_V(12, 12, 0.0);

        svd_decompose(matrix_P, matrix_U, matrix_S, matrix_V);

        //std::cout << "Matrix U (SVD decomposition) \n" << matrix_U << std::endl;
        //std::cout << "Matrix S (SVD decomposition) \n" << matrix_S << std::endl;
        //std::cout << "Matrix V (SVD decomposition) \n" << matrix_V << std::endl;

        // Now let's check if the SVD result is correct
        // Check 1: U is orthogonal, so U * U^T must be identity
        std::cout << "SVD Check 1 \n" << "U*U^T: \n" << matrix_U * transpose(matrix_U) << std::endl;
        // Check 2: V is orthogonal, so V * V^T must be identity
        std::cout << "SVD Check 2 \n" << "V*V^T: \n" << matrix_V * transpose(matrix_V) << std::endl;
        // Check 3: S must be a diagonal matrix
        std::cout << "SVD Check 3 \n" << "S: \n" << matrix_S << std::endl;
        // Check 4: according to the definition, A = U * S * V^T
        std::cout << "SVD Check 4 \n" << "M - U * S * V^T: \n" << matrix_P - matrix_U * matrix_S * transpose(matrix_V) << std::endl;

        //Computing vector m
        Vector vector_m = matrix_V.get_column(matrix_V.cols()-1);
        std::cout << "Vector m \n" << vector_m << std::endl;

        //Computing matrix M
        Matrix matrix_M(3, 4, 0.0);
        matrix_M.set_row(0, {vector_m[0], vector_m[1], vector_m[2], vector_m[3]});
        matrix_M.set_row(1, {vector_m[4], vector_m[5], vector_m[6], vector_m[7]});
        matrix_M.set_row(2, {vector_m[8], vector_m[9], vector_m[10], vector_m[11]});

        std::cout << "\nMatrix M\n" << matrix_M << std::endl;

        //Extracting intrinsic and extrinsic parameters from matrix M
        //Computing matrix A and vector b from matrix M
        Matrix matrix_A(3, 3, 0.0);
        Matrix vector_b(3, 1, 0.0);

        for (int i = 0; i < 3; i++) {
            matrix_A.set_column(i, matrix_M.get_column(i));
        }

        vector_b.set_column(0, matrix_M.get_column(3));

        std::cout << "Matrix A\n" << matrix_A << std::endl;
        std::cout << "Vector b\n" << vector_b << std::endl;

        //Computing intrinsic parameters
        Vector a1 = matrix_A.get_row(0);
        Vector a2 = matrix_A.get_row(1);
        Vector a3 = matrix_A.get_row(2);

        double rho = 1/norm(a3);
        double u0 = pow(rho, 2.0) * dot(a1, a3);
        double v0 = pow(rho, 2.0) * dot(a2, a3);
        double theta = acos(-1*dot(cross(a1, a3), cross(a2, a3))/(norm(cross(a1, a3))*norm(cross(a2, a3))));
        double alpha = pow(rho, 2.0)*norm(cross(a1, a3))*sin(theta);
        double beta = pow(rho, 2.0)*norm(cross(a2, a3))*sin(theta);

        std::cout << "Intrinsic parameters" << std::endl;
        std::cout << "rho = " << rho << std::endl;
        std::cout << "u0 = " << u0 << std::endl;
        std::cout << "v0 = " << v0 << std::endl;
        std::cout << "theta = " << theta << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "beta = " << beta << std::endl;

        //Computing extrinsic parameters
        Vector r1 = cross(a2, a3)/norm(cross(a2, a3));
        Vector r3 = rho * a3;
        Vector r2 = cross(r3,r1);

        Matrix K (3, 3, 0.0);
        K.set_row(0, {alpha, -alpha*(cos(theta)/sin(theta)), u0});
        K.set_row(1, {0, beta/sin(theta), v0});
        K.set_row(2, {0, 0, 1});
        Matrix invK;
        inverse(K, invK);
        Matrix small_t = rho*invK*vector_b;

        std::cout << "Extrinsic parameters" << std::endl;
        std::cout << "r1 = " << r1 << std::endl;
        std::cout << "r2 = " << r2 << std::endl;
        std::cout << "r3 = " << r3 << std::endl;
        std::cout << "K\n" << K << std::endl;
        std::cout << "small_t\n" << small_t << std::endl;

        //Computing output parameters
        fx = K[0][0];
        fy = K[1][1];
        cx = K[0][2];
        cy = K[1][2];
        skew = K[0][1];

        Matrix r_prev (3, 3, 0.0);
        r_prev.set_row(0, r1);
        r_prev.set_row(1, r2);
        r_prev.set_row(2, r3);
        R = Matrix33(r_prev);
        t = Vector3D(small_t[0][0], small_t[1][0], small_t[2][0]);

        std::cout << "Output parameters" << std::endl;
        std::cout << "fx = " << fx << std::endl;
        std::cout << "fy = " << fy << std::endl;
        std::cout << "cx = " << cx << std::endl;
        std::cout << "cy = " << cy << std::endl;
        std::cout << "skew = " << skew << std::endl;
        std::cout << "R\n" << R << std::endl;
        std::cout << "t\n" << t << std::endl;

        //Accuracy
        Matrix RT = Matrix(3, 4, 0.0);
        RT.set_row(0, {R[0][0], R[0][1], R[0][2], small_t[0][0]});
        RT.set_row(1, {R[1][0], R[1][1], R[1][2], small_t[1][0]});
        RT.set_row(2, {R[2][0], R[2][1], R[2][2], small_t[2][0]});

        Matrix M_check = K*RT;

        for (int i = 0; i < points_3d.size(); i++) {
            Matrix points3D;
            points3D = Matrix(4, 1, 0.0);
            points3D.set_column(0, {points_3d[i][0], points_3d[i][1], points_3d[i][2], 1});
            Matrix check_2d = M_check*points3D;

            std::cout << "\nAccuracy point " << i+1 << std::endl;
            std::cout << "3D coordinates\n" << points3D;
            std::cout << "Computed 2D coordinates\n" << check_2d/check_2d[2][0] << std::endl;
        }

        std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                     "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                     "\t\tif your calibration is successful or not.\n\n" << std::flush;

    } else {
        std::cout << "Input is invalid (number of correspondences < 6 or sizes of 2D/3D points don't match" << std::endl;

    }

    return true;
}