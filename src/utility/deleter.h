/*! \file deleter.h
    \brief gsl_interp_accelとgsl_splineのデリータを宣言・定義したヘッダファイル

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/


#ifndef _DELETER_H_
#define _DELETER_H_

#pragma once

#include <boost/checked_delete.hpp>     // for boost::checked_delete
#include <gsl/gsl_integration.h>    // for gsl_integration_glfixed_table_free
#include <xc.h>                         // for xc_func_end

namespace utility {
    //! A lambda expression.
    /*!
        gsl_integration_glfixed_tableへのポインタを解放するラムダ式
        \param ptable gsl_integration_glfixed_tableへのポインタ
    */
    static auto const gsl_integration_glfixed_table_deleter = [](auto ptable) noexcept {
        gsl_integration_glfixed_table_free(ptable);
    };

    //! A lambda expression.
    /*!
        xc_func_typeへのポインタを解放するラムダ式
        \param xcfunc xc_func_type へのポインタ
    */
    static auto const xcfunc_deleter = [](auto * xcfunc) {
        xc_func_end(xcfunc);
        boost::checked_delete(xcfunc);
    };
}

#endif  // _DELETER_H_
