#include "keys.h"
#include "fhe/common/mod_arith.h"
#include "fhe/common/ntt.h"
#include "fhe/common/rns_transform.h"

namespace hehub {

RlweKsk::RlweKsk(const RlweSk &sk_curr, const RlweSk &sk_orig,
                 const u64 additional_mod) {
    auto sk_curr_extended = sk_curr;
    sk_curr_extended.add_components({additional_mod});
    // no need to set 0 manually unless debugging mem since will be mult'd by p
#ifdef HEHUB_DEBUG
    std::fill(sk_curr_extended.last()->begin(), sk_curr_extended.last()->end(),
              (u64)0);
#endif

    // To encapsulate
    auto extended_moduli = sk_orig.modulus_vec();
    extended_moduli.push_back(additional_mod);
    auto sk_orig_extended(sk_orig);
    intt_negacyclic_inplace_lazy(sk_orig_extended);
    auto extended_part = rns_base_transform(sk_orig_extended, {additional_mod});
    sk_orig_extended.add_components({additional_mod});
    *sk_orig_extended.last() = std::move(extended_part[0]);
    ntt_negacyclic_inplace_lazy(sk_orig_extended);

    auto orig_components = sk_orig.component_count();
    std::vector<std::vector<u64>> rns_composition_basis(orig_components);
    for (size_t i = 0; i < orig_components; i++) {
        rns_composition_basis[i].resize(orig_components + 1, 0);
        rns_composition_basis[i][i] = additional_mod % extended_moduli[i];
    }
    *this = RgswCt(rgsw_encrypt_montgomery(sk_curr_extended, sk_orig_extended, 
                                           rns_composition_basis));
}

} // namespace hehub
