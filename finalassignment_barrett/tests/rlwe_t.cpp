#include "catch2/catch.hpp"
#include "fhe/primitives/rlwe.h"
#include <numeric>

using namespace hehub;

TEST_CASE("rlwe") {
    RnsPolyParams params = create_params(4096, {30, 30});
    RlwePt pt(params);
    RlweSk sk(params);

    // plaintext data
    const i64 DATUM_TEST = 123456;
    for (auto &component_poly : pt) {
        for (auto &datum : component_poly) {
            datum = DATUM_TEST;
        }
    }

    // encrypt & decrypt
    RlweCt ct = encrypt_core(pt, sk);
    RlwePt pt_recovered = decrypt_core(ct, sk);

    // check
    auto component_count = pt.component_count();
    REQUIRE(pt_recovered.component_count() == component_count);
    auto check_if_close = [=](auto status, auto coeff) {
        return status && (std::abs(i64(coeff) - DATUM_TEST) < 20);
    };
    for (auto &component_poly : pt_recovered) {
        REQUIRE(std::accumulate(component_poly.begin(), component_poly.end(),
                                true, check_if_close));
    }
}
