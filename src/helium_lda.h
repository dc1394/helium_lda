/*! \file helium_lda.h
    \brief VWN-LDAを用い、Kohn-Sham法でヘリウム原子のエネルギーを計算するクラスの宣言
    Copyright © 2019 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _HELIUM_LDA_H_
#define _HELIUM_LDA_H_

#pragma once

#include "utility/deleter.h"
#include <cstdint>                  // for std::int32_t
#include <memory>                   // for std::shared_ptr, std::unique_ptr            
#include <optional>                 // for std::optional
#include <valarray>                 // for std::valarray
#include <boost/multi_array.hpp>    // for boost::multi_array
#include <Eigen/Core>               // for Eigen::MatrixXd, Eigen::VectorXd

namespace helium_lda {
    //! A class.
    /*!
        VWN-LDAを用い、Kohn-Sham法でヘリウム原子のエネルギーを計算するクラス
    */
    class Helium_LDA final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
        */
        Helium_LDA();

        //! A destructor.
        /*!
            デストラクタ
        */
        ~Helium_LDA() = default;

        // #region publicメンバ関数

        //! A public member function.
        /*!
            SCF計算を行う
            \return SCF計算が正常に終了した場合はエネルギーを、しなかった場合はstd::nulloptを返す
        */
        std::optional<double> do_scfloop();

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            nalpha個のGTOによるヘリウム原子のエネルギーを計算する
            \param ep 一般化固有値問題のエネルギー固有値E'
            \return ヘリウム原子のエネルギー
        */
        double calc_energy(double ep);

        //! A private member function.
        /*!
            nalpha個のGTOによるヘリウム原子の交換相関エネルギーを計算する
            \return ヘリウム原子の交換相関エネルギー
        */
        double calc_exc_energy();

        //! A private member function.
        /*!
            K'計算用のgsl_functionを初期化する
        */
        void init_gsl_function_Kp();
        
        //! A private member function.
        /*!
            Kpq計算用のgsl_functionを初期化する
        */
        void init_gsl_function_Kpq();

        //! A private member function.
        /*!
            使用するGTOの数をユーザに入力させる
        */
        void input_nalpha();

        //! A private member function.
        /*!
            GTOの肩の係数が格納された配列を生成する
        */
        void make_alpha();

        //! A private member function.
        /*!
            全ての要素が、引数で指定された値で埋められたnalpha次元ベクトルを生成する
            \param val 要素を埋める値
        */
        void make_c(double val);

        //! A private member function.
        /*!
            交換相関積分が格納された、nalpha×nalphaの2次元配列を生成する
        */
        void make_exchcorrinteg();

        //! A private member function.
        /*!
            nalphaの数で、固有ベクトル、1電子積分および2電子積分からFock行列を生成する
        */
        void make_fockmatrix();

        //! A private member function.
        /*!
            1電子積分が格納された、nalpha×nalphaの2次元配列を生成する
        */
        void make_oneelectroninteg();

        //! A private member function.
        /*!
            nalpha次正方行列の重なり行列を生成する
        */
        void make_overlapmatrix();

        //! A private member function.
        /*!
            2電子積分が格納されたnalpha×nalpha×nalpha×nalphaの4次元配列を生成する
        */
        void make_twoelectroninteg();

        //! A private member function.
        /*!
            固有ベクトルを正規化する
        */
        void normalize();

        // #endregion privateメンバ関数

        // #region メンバ変数

        //! A private member variable (constant expression).
        /*!
            Gauss-Legendre積分の分点
        */
        static auto constexpr INTEGTABLENUM = 100;

        //! A private member variable (constant expression).
        /*!
            バッファサイズの上限
        */
        static auto constexpr MAXBUFSIZE = 32;

        //! A private member variable (constant expression).
        /*!
            SCF計算のループの上限
        */
        static auto constexpr MAXITER = 1000;

        //! A private member variable (constant expression).
        /*!
            積分区間の上限
        */
        static auto constexpr MAXR = 10.0;

        //! A private member variable (constant expression).
        /*!
            SCF計算のループから抜ける際のエネルギーの差の閾値
        */
        static auto constexpr SCFTHRESHOLD = 1.0E-15;

        //! A private member variable.
        /*!
            GTOの肩の係数が格納されたstd::vector
        */
        std::valarray<double> alpha_;

        //! A private member variable.
        /*!
            固有ベクトルC
        */
        Eigen::VectorXd c_;

        //! A private member variable.
        /*!
            Fock行列
        */
        Eigen::MatrixXd f_;
        
        //! A private member variable.
        /*!
            gslによるK'計算用の積分オブジェクト
        */
        gsl_function func_calc_Kp_;

        //! A private member variable.
        /*!
            gslによるKpq計算用の積分オブジェクト
        */
        gsl_function func_calc_Kpq_;

        //! A private member variable.
        /*!
            1電子積分が格納された2次元配列
        */
        boost::multi_array<double, 2> h_;

        //! A private member variable.
        /*!
            交換相関積分が格納された4次元配列
        */
        boost::multi_array<double, 2> k_;

        //! A private member variable.
        /*!
            使用するGTOの数
        */
        std::int32_t nalpha_ = 0;

        //! A private member variable (constant).
        /*!
            相関汎関数へのスマートポインタ
        */
        std::shared_ptr<xc_func_type> const pcfunc_;

        //! A private member variable (constant).
        /*!
            gsl_integration_glfixed_tableへのスマートポインタ
        */
        std::unique_ptr<gsl_integration_glfixed_table, decltype(utility::gsl_integration_glfixed_table_deleter)> const ptable_;

        //! A private member variable (constant).
        /*!
            交換汎関数へのスマートポインタ
        */
        std::shared_ptr<xc_func_type> const pxfunc_;

        //! A private member variable.
        /*!
            2電子積分が格納された4次元配列
        */
        boost::multi_array<double, 4> q_;

        //! A private member variable.
        /*!
            重なり行列
        */
        Eigen::MatrixXd s_;

        //! A private member variable.
        /*!
            交換相関積分が格納された2次元配列
        */
        boost::multi_array<double, 2> xc_;

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A public copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        Helium_LDA(Helium_LDA const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Helium_LDA & operator=(Helium_LDA const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
        
}

#endif  // HELIUM_LDA_H_
