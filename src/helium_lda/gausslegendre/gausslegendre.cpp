/*! \file gausslegendre.cpp
    \brief Gauss-Legendre積分を行うクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "../../alglib/src/integration.h"	// for alglib
#include "gausslegendre.h"
#include <stdexcept>					    // for std::runtime_error
#include <boost/cast.hpp>				    // for boost::numeric_cast

namespace gausslegendre {
    // #region コンストラクタ

    Gauss_Legendre::Gauss_Legendre(std::int32_t n) : n_(n)
    {
        alglib::ae_int_t info = 0;
        alglib::real_1d_array x, w;

        alglib::gqgenerategausslegendre(
            boost::numeric_cast<alglib::ae_int_t>(n),
            info,
            x,
            w);

        if (info != 1) {
            throw std::runtime_error("alglib::gqgenerategausslegendreが失敗");
        }

        x_.assign(x.getcontent(), x.getcontent() + x.length());
        w_.assign(w.getcontent(), w.getcontent() + w.length());
    }

    // #endregion コンストラクタ
}

