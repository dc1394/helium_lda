/*! \file helium_lda.cpp
    \brief VWN-LDAを用い、Kohn-Sham法でヘリウム原子のエネルギーを計算するクラスの実装
    Copyright © 2019 @dc1394 All Rights Reserved.
    (but this is originally adapted by Paolo Giannozzi for helium_hf_gauss.c from http://www.fisica.uniud.it/~giannozz/Corsi/MQ/Software/C/helium_hf_gauss.c )
    This software is released under the BSD 2-Clause License.
*/

#include "helium_lda.h"
#include "utility/deleter.h"
#include <cmath>                                // for std::pow, std::sqrt
#include <iostream>                             // for std::cerr, std::cin, std::cout
#include <memory>                               // for std::unique_ptr
#include <optional>                             // for std::make_optional, std::nullopt
#include <utility>                              // for std::make_pair
#include <boost/assert.hpp>                     // for BOOST_ASSERT
#include <boost/format.hpp>                     // for boost::format
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi
#include <Eigen/Eigenvalues>                    // for Eigen::GeneralizedSelfAdjointEigenSolver
#include <gsl/gsl_integration.h>

namespace helium_lda {
    Helium_LDA::Helium_LDA()
        :   pcfunc_(new xc_func_type, utility::xcfunc_deleter),
            pxfunc_(new xc_func_type, utility::xcfunc_deleter)
    {
        xc_func_init(pcfunc_.get(), XC_LDA_C_VWN, XC_POLARIZED);
        xc_func_init(pxfunc_.get(), XC_LDA_X, XC_POLARIZED);

        // 使用するGTOの数を入力
        input_nalpha();

        f_ = Eigen::MatrixXd::Zero(nalpha_, nalpha_);
        h_.resize(boost::extents[nalpha_][nalpha_]);
        q_.resize(boost::extents[nalpha_][nalpha_][nalpha_][nalpha_]);
        s_ = Eigen::MatrixXd::Zero(nalpha_, nalpha_);
    }

    std::optional<double> Helium_LDA::do_scfloop()
    {
        // GTOの肩の係数が格納された配列を生成
        make_alpha();

        // 1電子積分が格納された2次元配列を生成
        make_oneelectroninteg();

        // 2電子積分が格納された4次元配列を生成
        make_twoelectroninteg();

        // 重なり行列を生成
        make_overlapmatrix();

        // 全て1.0で初期化された固有ベクトルを生成
        make_c(1.0);

        // 固有ベクトルを正規化
        normalize();

        // 新しく計算されたエネルギー
        auto enew = 0.0;

        // SCFループ
        for (auto iter = 1; iter < Helium_LDA::MAXITER; iter++) {
            // Fock行列を生成
            make_fockmatrix();

            // 一般化固有値問題を解く
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(f_, s_);

            // E'を取得
            auto const ep = es.eigenvalues()[0];

            // 固有ベクトルを取得
            c_ = es.eigenvectors().col(0);

            // 固有ベクトルを正規化
            normalize();

            // 前回のSCF計算のエネルギーを保管
            auto const eold = enew;

            // 今回のSCF計算のエネルギーを計算する
            enew = getenergy(ep);

            std::cout << boost::format("Iteration # %2d: HF eigenvalue = %.14f, energy = %.14f\n") % iter % ep % enew;

            // SCF計算が収束したかどうか
            if (std::fabs(enew - eold) < Helium_LDA::SCFTHRESHOLD) {
                // 収束したのでそのエネルギーを返す
                return std::make_optional(enew);
            }
        }

        // SCF計算が収束しなかった
        return std::nullopt;
    }

    double Helium_LDA::getenergy(double ep) const
    {
        auto e = ep;
        for (auto p = 0; p < nalpha_; p++) {
            for (auto q = 0; q < nalpha_; q++) {
                // E = E' + Cp * Cq * hpq
                e += c_[p] * c_[q] * h_[p][q];
            }
        }

        return e;
    }

    void Helium_LDA::input_nalpha()
    {
        while (true) {
            std::cout << "使用するGTOの個数を入力してください (3, 4 or 6): ";
            std::cin >> nalpha_;

            if (!std::cin.fail() && (nalpha_ == 3 || nalpha_ == 4 || nalpha_ == 6)) {
                break;
            }

            std::cin.clear();
            std::cin.ignore(Helium_LDA::MAXBUFSIZE, '\n');
        }
    }

    void Helium_LDA::make_alpha()
    {
        switch (nalpha_) {
        case 3:
            alpha_ = { 0.31364978999999998, 1.1589229999999999, 6.3624213899999997 };
            break;

        case 4:
            alpha_ = { 0.297104, 1.236745, 5.749982, 38.2166777 };
            break;

        case 6:
            alpha_ = { 0.100112428, 0.24307674700000001, 0.62595526599999995, 1.8221429039999999, 6.5131437249999999, 35.523221220000003 };
            break;

        default:
            BOOST_ASSERT(!"switch文のdefaultに来てしまった！");
            break;
        }
    }

    void Helium_LDA::make_c(double val)
    {
        c_.resize(nalpha_);
        // 固有ベクトルCの要素を全てvalで初期化
        for (auto i = 0; i < nalpha_; i++) {
            c_[i] = val;
        }
    }

    void Helium_LDA::make_fockmatrix()
    {
        for (auto p = 0; p < nalpha_; p++) {
            for (auto qi = 0; qi < nalpha_; qi++) {
                // Fpq = hpq + ΣCr * Cs * Qprqs
                f_(p, qi) = h_[p][qi];

                for (auto r = 0; r < nalpha_; r++) {
                    for (auto s = 0; s < nalpha_; s++) {
                        f_(p, qi) += c_[r] * c_[s] * q_[p][r][qi][s];
                    }
                }
            }
        }
    }

    void Helium_LDA::make_oneelectroninteg()
    {
        using namespace boost::math::constants;

        for (auto p = 0; p < nalpha_; p++) {
            for (auto q = 0; q < nalpha_; q++) {
                // αp + αq
                auto const appaq = alpha_[p] + alpha_[q];

                // hpq = 3αpαqπ^1.5 / (αp + αq)^2.5 - 4π / (αp + αq)
                h_[p][q] = 3.0 * alpha_[p] * alpha_[q] * std::pow((pi<double>() / appaq), 1.5) / appaq -
                    4.0 * pi<double>() / appaq;
            }
        }
    }

    void Helium_LDA::make_overlapmatrix()
    {
        using namespace boost::math::constants;

        for (auto p = 0; p < nalpha_; p++) {
            for (auto q = 0; q < nalpha_; q++) {
                // Spq = (π / (αp + αq))^1.5
                s_(p, q) = std::pow((pi<double>() / (alpha_[p] + alpha_[q])), 1.5);
            }
        }
    }

    void Helium_LDA::make_twoelectroninteg()
    {
        using namespace boost::math::constants;

        for (auto p = 0; p < nalpha_; p++) {
            for (auto qi = 0; qi < nalpha_; qi++) {
                for (auto r = 0; r < nalpha_; r++) {
                    for (auto s = 0; s < nalpha_; s++) {
                        // Qprqs = 2π^2.5 / [(αp + αq)(αr + αs)√(αp + αq + αr + αs)]
                        q_[p][r][qi][s] = 2.0 * std::pow(pi<double>(), 2.5) /
                            ((alpha_[p] + alpha_[qi]) * (alpha_[r] + alpha_[s]) *
                                std::sqrt(alpha_[p] + alpha_[qi] + alpha_[r] + alpha_[s]));
                    }
                }
            }
        }
    }

    void Helium_LDA::normalize()
    {
        using namespace boost::math::constants;

        std::unique_ptr<gsl_integration_glfixed_table, decltype(utility::gsl_integration_glfixed_table_deleter)> const
            ptable(gsl_integration_glfixed_table_alloc(Helium_LDA::INTEGTABLENUM), utility::gsl_integration_glfixed_table_deleter);
        gsl_function F;

        F.function = [](double x, void * params)
        {
            auto const [alpha, c] = *(reinterpret_cast<std::pair< std::valarray<double>, Eigen::VectorXd> *>(params));
            auto const nalpha = alpha.size();

            auto f = 0.0;
            for (auto p = 0U; p < nalpha; p++)
            {
                f += c[p] * std::exp(-alpha[p] * x * x);
            }

            return 4.0 * pi<double>() * x * x * f * f;
        };

        auto params = std::make_pair(alpha_, c_);
        F.params = reinterpret_cast<void *>(&params);

        auto const sum = gsl_integration_glfixed(&F, 0.0, 10.0, ptable.get());

        for (auto p = 0; p < nalpha_; p++)
        {
            c_[p] /= sum;
        }
    }
}
