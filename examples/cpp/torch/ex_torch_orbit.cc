/**
 * @file ex_torch_orbit.cc
 * @author Stanford NAV LAB
 * @brief  Compute the state transition matrix (STM) using Pytorch autograd (as potential alternative to autodiff)
 * @version 0.1
 * @date 2024-12-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <lupnt/lupnt.h>
#include <chrono>
#include <cmath>
#include <torch/torch.h>

using namespace torch::autograd;
using namespace lupnt;

using ODET = std::function<torch::Tensor(double t, torch::Tensor& x)>;

/******** Autodiff Functions *******************/
VecX crtbp_ad(Real t, const VecX& x){
    double mu = 0.0121505;
    Real x1 = x(0);
    Real x2 = x(1);
    Real x3 = x(2);
    Real x4 = x(3);
    Real x5 = x(4);
    Real x6 = x(5);
    
    Real r1 = sqrt(pow(x1 + mu, 2) + pow(x2, 2) + pow(x3, 2));
    Real r2 = sqrt(pow(x1 - 1 + mu, 2) + pow(x2, 2) + pow(x3, 2));
    
    Real dx1 = x4;
    Real dx2 = x5;
    Real dx3 = x6;
    Real dx4 = 2 * x5 + x1 - (1 - mu) * (x1 + mu) / pow(r1, 3) - mu * (x1 - 1 + mu) / pow(r2, 3);
    Real dx5 = -2 * x4 + x2 - (1 - mu) * x2 / pow(r1, 3) - mu * x2 / pow(r2, 3);
    Real dx6 = -(1 - mu) * x3 / pow(r1, 3) - mu * x3 / pow(r2, 3);
    
    VecX dx(6);
    dx << dx1, dx2, dx3, dx4, dx5, dx6;
    
    return dx;
}

 
/******** Torch Functions *******************/
torch::Tensor crtbp_torch(double t, const torch::Tensor& x) {
  auto mu = 0.0121505;
  auto x1 = x[0];
  auto x2 = x[1];
  auto x3 = x[2];
  auto x4 = x[3];
  auto x5 = x[4];
  auto x6 = x[5];

  auto r1 = torch::sqrt(torch::pow(x1 + mu, 2) + torch::pow(x2, 2) + torch::pow(x3, 2));
  auto r2 = torch::sqrt(torch::pow(x1 - 1 + mu, 2) + torch::pow(x2, 2) + torch::pow(x3, 2));

  auto dx1 = x4;
  auto dx2 = x5;
  auto dx3 = x6;
  auto dx4 = 2 * x5 + x1 - (1 - mu) * (x1 + mu) / torch::pow(r1, 3) - mu * (x1 - 1 + mu) / torch::pow(r2, 3);
  auto dx5 = -2 * x4 + x2 - (1 - mu) * x2 / torch::pow(r1, 3) - mu * x2 / torch::pow(r2, 3);
  auto dx6 = -(1 - mu) * x3 / torch::pow(r1, 3) - mu * x3 / torch::pow(r2, 3);

  torch::Tensor dx = torch::stack({dx1, dx2, dx3, dx4, dx5, dx6});

  return dx;
}

torch::Tensor rk4_torch(const ODET& f, double t, torch::Tensor& x, double dt) {
  // Evaluate `f` (i.e., `dx`) at the 4 locations defined by the RK4 method
  auto k_1 = f(t, x) * dt;

  auto t1 = t + dt / 2.0;
  auto x1 = x + k_1 / 2.0;
  auto k_2 = f(t1, x1) * dt;

  auto t2 = t + dt / 2.0;
  auto x2 = x + k_2 / 2.0;
  auto k_3 = f(t2, x2) * dt;

  auto t3 = t + dt;
  auto x3 = x + k_3;
  auto k_4 = f(t3, x3) * dt;

  // Average the 4 derivatives to approximate `dx`
  auto dx = (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) / 6.0;

  return x + dx;
}

struct STMResultT {
    torch::Tensor final_state;  // shape [6]
    torch::Tensor stm;          // shape [6, 6]
};

STMResultT propagate_with_stm_torch(torch::Tensor& x_init,
                            double t0,
                            double dt,
                            int num_steps)
{
    // 1) Clone and enable gradient for initial state:
    // Make sure x_init is float or double (not int).
    x_init.set_requires_grad(true);
    auto x = x_init.clone();

    // 2) Integrate from t0 to t0 + num_steps*dt
    double time = t0;
    for (int step = 0; step < num_steps; ++step) {
        x = rk4_torch(crtbp_torch, time, x, dt);
        time += dt;
    }

    // 3) Build the 6x6 STM by doing component-wise backward
    const int n = x.size(0); // should be 6
    torch::Tensor stm = torch::zeros({n, n}, x.options());

    for (int i = 0; i < n; ++i) {
        // Zero grads on x_init each time
        // We'll read x_init.grad() after backward
        // to get partial derivatives of x_final[i] wrt x_init
        if (x_init.grad().defined()) {
            x_init.grad().zero_();
        }

        // If not the last iteration, we do retain_graph(true)
        // so the graph stays alive for subsequent calls
        bool retain = (i < (n-1));

        // Now backprop from x_final[i]
        // This tells autograd: "d(x_final[i]) / d(x_init)"
        // We pass a ones() gradient of shape [], meaning "1.0" for that scalar.
        x[i].backward({}, retain);

        // Now x_init.grad() is a [6] tensor holding the partial derivatives
        // of x_final[i] wrt x_init[j].
        auto grad_i = x_init.grad();

        // Place this row into STM:
        // row i of STM is the gradient of x_final[i] w.r.t x_init
        // shape: [6]
        stm[i] = grad_i;
    }

    STMResultT result;
    result.final_state = x;
    result.stm         = stm;
    return result;
}

int main() {
    // Suppose a sample initial state:
    // x_init = [x1, x2, x3, x4, x5, x6]
    Real t0 = 0.0;
    Real dt = 0.01;
    int num_steps = 100;
    Real tf = t0 + num_steps * dt;

    // Compute final state & STM using PyTorch autograd -------------------------
    std::cout << "--------------------------------------------------\n";
    // clock start
    torch::Tensor x_init = torch::tensor({1.2, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        torch::TensorOptions().dtype(torch::kFloat64));

    auto start = std::chrono::high_resolution_clock::now();
    STMResultT res = propagate_with_stm_torch(x_init, t0.val(), dt.val(), num_steps);

    std::cout << "Final state after " << num_steps << " steps:  (Torch)\n"
              << res.final_state << "\n\n";

    std::cout << "STM (dXfinal/dXinit):  (Torch) \n" << res.stm << "\n";

    // clock end
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Elapsed time (Pytorch): " << elapsed.count() << " s\n\n";

    std::cout << "--------------------------------------------------\n\n";


    // Compute final state & STM using Autodiff -------------------------

    // clock start
    VecX x_init_ad(6);
    x_init_ad << 1.2, 0.0, 0.0, 0.0, 0.0, 0.0;

    auto start_ad = std::chrono::high_resolution_clock::now();
    Ptr<NumericalPropagator> propagator = MakePtr<NumericalPropagator>(IntegratorType::RK4);

    MatXd stm_ad(6, 6);
    VecX x_final_ad = propagator->Propagate(crtbp_ad, t0, tf, x_init_ad, dt, &stm_ad);

    std::cout << "Final state after " << num_steps << " steps:  (Autodiff)\n"
              << x_final_ad.transpose() << "\n\n";
    
    std::cout << "STM (dXfinal/dXinit):  (Autodiff) \n" << stm_ad << "\n";

    // clock end
    auto end_ad = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ad = end_ad - start_ad;

    std::cout << "Elapsed time (Autodiff): " << elapsed_ad.count() << " s\n\n";

    std::cout << "--------------------------------------------------\n\n";

    return 0;
}