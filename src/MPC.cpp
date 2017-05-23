#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 8;
double dt = 0.15;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// reference values for errors and velocity
double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 100;

// Solver takes all state and actuator variables in a singular vector.
// Set up indeces for easier access later
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t deltaPsi_start = epsi_start + N;
size_t a_start = deltaPsi_start + N;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;
      
    // increment cost for each timestep's cross track and orientation error
    for (int i = 0; i < N; i++) {
        fg[0] += CppAD::pow(vars[cte_start+i] - ref_cte, 2);
        fg[0] += CppAD::pow(vars[epsi_start+1] - ref_epsi, 2);
        fg[0] += CppAD::pow(vars[v_start+1] - ref_v, 2);
    }
    // penalize big value gaps in sequential actuations
    for (int i = 0; i < N-2; i++) {
        fg[0] += CppAD::pow(vars[deltaPsi_start+i+1] - vars[deltaPsi_start+i], 2);
        fg[0] += CppAD::pow(vars[a_start+i+1] - vars[a_start+i], 2);
    }
    // penalize the use of actuators
    for (int i = 0; i < N-1; i++) {
        fg[0] += CppAD::pow(vars[deltaPsi_start+i], 2);
        fg[0] += CppAD::pow(vars[a_start+i], 2);
    }
    

    // Setup Constraints
    // We add 1 to each of the starting indices due to cost being located at index 0 of `fg`.
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[psi_start + 1] = vars[psi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> delta0 = vars[deltaPsi_start + i];
      AD<double> a0 = vars[a_start + i];
      AD<double> epsi0 = vars[epsi_start + i];
      
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i +1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];

      AD<double> f_x0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0,2));
      
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is so CppAD can compute derivatives and pass these to the solver.
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[2 + x_start + i]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y_start + i]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_start + i]  = psi1 - (psi0 + v0 / Lf * delta0 * dt);
      fg[2 + v_start + i]    = v1 - (v0 + a0 * dt);
      fg[2 + cte_start + i]  = cte1 - (f_x0 - y0 + v0 * CppAD::sin(epsi0) * dt);
      fg[2 + epsi_start + i] = epsi1 - (psi0 - psides0 + v0 * dt / Lf * dt);
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
  // element vector and there are 10 timesteps. The number of variables is: 4 * 10 + 2 * 9
  int dim_state_vec = 6; // x, y, psi, v, cte, epsi
  int dim_actuator_vec = 2; // steering, acceleration
  
  // number of independent variables
  // CHANGE... NOT N-1
  //size_t n_vars = N * dim_state_vec + (N - 1) * dim_actuator_vec;
  size_t n_vars = N * dim_state_vec + N * dim_actuator_vec;
  
  // TODO: Set the number of constraints
  size_t n_constraints = N * dim_state_vec;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  // set initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;
          
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  for (int i = 0; i < deltaPsi_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // max steering input: -25 to 25 degrees in radians
  for (int i = deltaPsi_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332; // try -0.2
    vars_upperbound[i] = 0.436332;
  }
  // max acceleration between -1 and 1
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // constrain initial state to input state
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

  
  
  // object that computes objective and constraints
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
 
  // save predicted x and y values
  x_vals.clear();
  y_vals.clear();
  for (int i=0; i<N; i++) {
    x_vals.push_back(solution.x[x_start+i]);
    y_vals.push_back(solution.x[y_start+i]);
  }
  
  for (int i=0; i<N; i++) {
      std::cout << "CTE: " << solution.x[cte_start+i] << " DELTA: " << solution.x[deltaPsi_start+i] << std::endl;
  }

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  return {solution.x[deltaPsi_start], solution.x[a_start]};
}