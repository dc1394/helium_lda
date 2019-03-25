/*! \file helium_lda_main.cpp
    \brief VWN-LDAを用い、Kohn-Sham法でヘリウム原子のエネルギーを計算する
    Copyright © 2019 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "helium_lda.h"
#include <iostream>             // for std::cerr, std::cout
#include <boost/format.hpp>     // for boost::format

int main()
{
    helium_lda::Helium_LDA hl;
    if (auto const res(hl.do_scfloop()); res) {
        std::cout << boost::format("SCF計算が収束しました: energy = %.14f (Hartree)") % (*res) << std::endl;

        return 0;
    }

    std::cerr << "SCF計算が収束しませんでした" << std::endl;
    return -1;
}
