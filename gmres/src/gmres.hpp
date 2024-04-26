#include "../../hessenberg/src/arnoldi.hpp"
#include "../../hessenberg/src/givens.hpp"

namespace slv {

    template<typename T>
    std::vector<T> gmres(const mtrx::csr<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T tolerance, int N_max) {
        std::vector<T> r0 = b - A*x0;
        std::vector<T> x = x0;
        T beta = std::sqrt(r0*r0);
        std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>> H_q = krylov::arnoldi_onb(A, r0, 1);
        std::vector<std::vector<T>> R = H_q.first; 
        R.reserve(N_max);
        T _a = R[0][0];
        T _b = R[0][1];
        T _r = std::sqrt(_a*_a + _b*_b);
        std::pair<T, T> c_s(_a/_r, -_b/_r);
        R[0][0] = _r;
        R[0][1] = 0;
        mtrx::dense<T> Omega(2, 2, {c_s.first, -c_s.second, c_s.second, c_s.first});
        T gg = beta*Omega(0,0);
        T yy = gg/R[0][0];
        x = x0 + yy*H_q.second[0];
        auto d = A*x - b;
        for (int i = 1; (i < N_max) && !(std::sqrt(d*d) < tolerance); i++) {
            krylov::arnoldi_iter(A, H_q.second, H_q.first);
            
            auto new_h = H_q.first.back();
            for (int i = 0; i < new_h.size()-1; i++) {
                new_h[i] = 0;
                for (int j = 0; j < Omega.width(); j++) {
                    new_h[i] += Omega(i, j)*H_q.first.back()[j];
                }
            }
            R.push_back(new_h);
            
            c_s = givens::zero_last(R, Omega);

            std::vector<T> new_Omega((Omega.height()+1)*(Omega.width()+1));
            for (int i = 0; i < Omega.height()-1; i++) {
                for (int j = 0; j < Omega.width()+1; j++) {
                    new_Omega[i*(Omega.width()+1)+j] = (j < Omega.width() ? Omega(i, j) : 0);
                }
            }
            i = Omega.height() - 1;
            for (int j = 0; j < Omega.width()+1; j++) {
                new_Omega[i*(Omega.width()+1)+j] = (j < Omega.width() ? c_s.first*Omega(i,j) : c_s.second);
            }
            i = Omega.height();
            for (int j = 0; j < Omega.width()+1; j++) {
                new_Omega[i*(Omega.width()+1)+j] = (j < Omega.width() ? -c_s.second*Omega(i-1,j) : c_s.first);
            }

            Omega = mtrx::dense(Omega.height()+1, Omega.width()+1, new_Omega);

            std::vector<T> g(R.size());
            for (int k = 0; k < g.size(); k++) {
                g[k] = beta*Omega(k, 0);
            }

            std::vector<T> y(R.size());        
            for (int k = y.size() - 1; k >= 0; k--) {
                T w = g[k];
                for (int j = y.size() - 1; j > k; j--) {
                    w -= y[j]*R[j][k];
                }
                y[k] = w/R[k][k];
            }

            std::vector<T> Q_y(H_q.second[0].size());
            for (int k = 0; k < Q_y.size(); k++) {
                T sum = 0;
                for (int j = 0; j < y.size(); j++) {
                    sum += y[j]*H_q.second[j][k];
                }
                Q_y[k] = sum;
            } 
            x = x0 + Q_y;
            d = A*x - b;
        }
        return x;
    }

}