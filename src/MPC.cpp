#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 25;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//

// Set vehicle speed.
const double ref_v = 60;

// Establish when one variable starts and another ends.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }
  
  AD<double> CalRcurv(Eigen::VectorXd coeffs, AD<double> x) {
    AD<double> dx_dy = 0.0;
    AD<double> d2x_dy2 = 0.0;
    AD<double> rcurv = 0.0;
    
    if (coeffs.size() > 2) {
      for (int i = 1; i < coeffs.size(); i++) {
        dx_dy += float(i) * coeffs[i] * CppAD::pow(x, i -1.0);
      }
    } else if (coeffs.size() == 2) {
      dx_dy = coeffs[1];
    }
    
    if (coeffs.size() > 3) {
      for (int i = 2; i < coeffs.size(); i++) {
        d2x_dy2 += float(i) * (float(i) - 1.0) * coeffs[i] * CppAD::pow(x, i - 2.0);
      }
    } else if (coeffs.size() == 3) {
      d2x_dy2 = coeffs[2];
    }
    
    if (d2x_dy2 != 0) {
      rcurv = CppAD::pow((1 + CppAD::pow(dx_dy, 2.0)), 3.0/2.0) / CppAD::abs(d2x_dy2);
    } else {
      rcurv = 1000.0;
    }
    
    return rcurv;
  }
  
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // fg a vector of constraints, x is a vector of constraints.
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    
    // Initialize cost.
    // The cost is stored is the first element of `fg`.
    fg[0] = 0;
    
    // Reference State Cost
    for (int i = 0; i < N; i++) {
      AD<double> rcurv;
      AD<double> d;
      
      rcurv = FG_eval::CalRcurv(coeffs, vars[x_start + i]);
      if (rcurv >= 1000) {
        d = 0.0;
      } else {
        d = 1.0 / rcurv;
      }
      AD<double> epsi = vars[epsi_start + i];
      AD<double> cte = vars[cte_start + i];
      AD<double> verr;
      if (10.0 * CppAD::abs(d) < ref_v) {
        verr = vars[v_start + i] - (ref_v - 0.0);
      } else {
        verr = vars[v_start + i];
      }
      
      fg[0] += 1.0 * CppAD::pow(epsi, 2);
      fg[0] += 1.0 * CppAD::pow(cte, 2);
      fg[0] += 0.1 * CppAD::pow(verr, 2);
    }
    
    // Define the cost related actuators.
    
    for (int i = 0; i < N -1; i++) {
        AD<double> delta = vars[delta_start + i];
        AD<double> a = vars[a_start + i];
        /*
        fg[0] += CppAD::pow(a, 2);
        */
        fg[0] += CppAD::pow(delta,2);
    }
    
    
    // Define the cost related actuators derivatives.
    for (int i = 0; i < N -2; i++) {
        AD<double> d_delta = vars[delta_start + i + 1] - vars[delta_start + i];
        AD<double> d_a = vars[a_start + i + 1] - vars[a_start + i];
        
        fg[0] += 500.0 * CppAD::pow(d_delta, 2);
        fg[0] += CppAD::pow(d_a, 2);
    }
    
    //
    // Setup Constraints
    //
    
    // Initial constraints
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    
    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];
      
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> delta0 = vars[delta_start + i];
      AD<double> a0 = vars[a_start + i];
      AD<double> f0;
      for (int i = 0; i < coeffs.size(); i++) {
        f0 += coeffs[i] *pow(x0, i);
      }
      AD<double> psides0;
      if (coeffs.size() > 2) {
        for (int i = 1; i < coeffs.size(); i++) {
          psides0 += coeffs[i] * pow(x0, i - 1) * i;
        }
      } else {
        psides0 = coeffs[1];
      }
      
      AD<double> d_psi0 = psi0 - atan(psides0);
      if (d_psi0 < -M_PI) {
        d_psi0 += 2 * M_PI;
      } else if (d_psi0 > M_PI) {
        d_psi0 -= 2 * M_PI;
      }
      
      fg[2 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_start + i] = psi1 -(psi0 + v0 / Lf * delta0 * dt);
      fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte_start + i] = cte1 - ((f0 - y0) + v0 * CppAD::sin(d_psi0) * dt);
      fg[2 + epsi_start + i] = epsi1 - (d_psi0 + v0 / Lf * delta0 * dt);
    }
    
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];
  
  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6 * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;
  
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  
  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -100.0;
    vars_upperbound[i] = 100.0;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  // limit initial constraint to initial state.
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);
  
  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);
  for (int i = 0; i < 10; i++) {
    result.push_back(solution.x[x_start + i]);
  }
  for (int i = 0; i < 10; i++) {
    result.push_back(solution.x[y_start + i]);
  }
  
  // 100msec latency state value.
  result.push_back(solution.x[psi_start + 2]);
  result.push_back(solution.x[v_start + 2]);
  result.push_back(solution.x[cte_start + 2]);
  result.push_back(solution.x[epsi_start + 2]);
  
  return result;
}
