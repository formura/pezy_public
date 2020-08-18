#pragma once

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#include <CL/cl2.hpp>
#include <CL/cl_pz_ext.h>

static auto clExtGetProfilePEStatistics = (pfnPezyExtGetProfilePEStatistics)clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_get_profile_pe_statistics");
static auto clExtGetProfilePE = (pfnPezyExtGetProfilePE)clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_get_profile_pe");
static auto clExtGetProfileCacheStatistics = (pfnPezyExtGetProfileCacheStatistics)clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_get_profile_cache_statistics");
static auto clExtGetProfileCache = (pfnPezyExtGetProfileCache)clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_get_profile_cache");

static auto clExtSetProfile = (pfnPezyExtSetProfile)clGetExtensionFunctionAddressForPlatform(nullptr, "pezy_set_profile");


#define PREF_MSTRD 0x40
#define PREF_MSTRDHIT 0x44
#define PREF_MSTWR 0x48
#define PREF_MSTWRHIT 0x4c

pzcl_profile_pe &operator-=(pzcl_profile_pe &l, const pzcl_profile_pe &r)
{
    l.run -= r.run;
    l.stall -= r.stall;
    l.wait -= r.wait;
    return l;
}

pzcl_profile_cache &operator-=(pzcl_profile_cache &l, const pzcl_profile_cache &r)
{
    l.read_request -= r.read_request;
    l.read_hit -= r.read_hit;
    l.write_request -= r.write_request;
    l.write_hit -= r.write_hit;
    return l;
}

class Profiles {
public:
    Profiles(cl::Context &ctx)
    {
        pe_stats = { sizeof(pzcl_profile_pe_stats) };
        clExtGetProfilePEStatistics(ctx(), 0, &pe_stats);

        for (size_t i = 0; i < pe_prof.size(); i++) {
            clExtGetProfilePE(ctx(), 0, i, &pe_prof[i]);
        }

        cache_stats_l1 = { sizeof(pzcl_profile_cache_stats) };
        cache_stats_l2 = { sizeof(pzcl_profile_cache_stats) };
        clExtGetProfileCacheStatistics(ctx(), 0, PZCL_EXT_PROFILE_CACHE_L1, &cache_stats_l1);
        clExtGetProfileCacheStatistics(ctx(), 0, PZCL_EXT_PROFILE_CACHE_L2, &cache_stats_l2);

        for (size_t i = 0; i < cache_prof_l1.size(); i++) {
            clExtGetProfileCache(ctx(), 0, PZCL_EXT_PROFILE_CACHE_L1, i, &cache_prof_l1[i]);
        }
        for (size_t i = 0; i < cache_prof_l2.size(); i++) {
            clExtGetProfileCache(ctx(), 0, PZCL_EXT_PROFILE_CACHE_L2, i, &cache_prof_l2[i]);
        }
    }

    Profiles &operator-=(const Profiles &r)
    {
        for (size_t i = 0; i < pe_prof.size(); i++) {
            pe_prof[i] -= r.pe_prof[i];
        }
        for (size_t i = 0; i < cache_prof_l1.size(); i++) {
            cache_prof_l1[i] -= r.cache_prof_l1[i];
        }
        for (size_t i = 0; i < cache_prof_l2.size(); i++) {
            cache_prof_l2[i] -= r.cache_prof_l2[i];
        }

        pe_stats.elapse_ns -= r.pe_stats.elapse_ns;

        return *this;
    }

    void print() const
    {
        using namespace std;

#if 1
        printf("PE profiles:\n");
        for (int i = 0; i < 2048; i++) {
            printf("  - { id: %d, run: %ld, stall: %ld(%02.2f%%), wait: %ld(%02.2f%%) }\n", i,
                pe_prof[i].run,
                pe_prof[i].stall,
                pe_prof[i].stall * 100.0 / pe_prof[i].run,
                pe_prof[i].wait,
                pe_prof[i].wait * 100.0 / pe_prof[i].run);
        }

        printf("\n");
        printf("Cache profiles:\n");
        printf("  L1:\n");
        for (int i = 0; i < (int)cache_prof_l1.size(); i++) {
            printf("    - { id: %d, read_req: %d, read_hit: %d, read_hit_ratio: %02.2f%%, write_req: %d, write_hit: %d, write_hit_ratio: %02.f%% }\n", i,
                cache_prof_l1[i].read_request,
                cache_prof_l1[i].read_hit,
                (double)cache_prof_l1[i].read_hit / cache_prof_l1[i].read_request * 100,
                cache_prof_l1[i].write_request,
                cache_prof_l1[i].write_hit,
                (double)cache_prof_l1[i].write_hit / cache_prof_l1[i].write_request * 100);
        }
        printf("  L2:\n");
        for (int i = 0; i < (int)cache_prof_l2.size(); i++) {
            printf("    - { id: %d, read_req: %d, read_hit: %d, read_hit_ratio: %02.2f%%, write_req: %d, write_hit: %d, write_hit_ratio: %02.f%% }\n", i,
                cache_prof_l2[i].read_request,
                cache_prof_l2[i].read_hit,
                (double)cache_prof_l2[i].read_hit / cache_prof_l2[i].read_request * 100,
                cache_prof_l2[i].write_request,
                cache_prof_l2[i].write_hit,
                (double)cache_prof_l2[i].write_hit / cache_prof_l2[i].write_request * 100);
        }
#endif

        uint64_t runs = 0, stalls = 0, waits = 0;
        for (auto &p : pe_prof) {
            runs += p.run;
            stalls += p.stall;
            waits += p.wait;
        }
        printf("\n");
        printf("PE avg statistics:\n");
        printf("  run:   %10ld\n", runs / pe_prof.size());
        printf("  stall: %10ld (%5.2f%%)\n", stalls / pe_prof.size(), stalls * 100.0 / runs);
        printf("  wait:  %10ld (%5.2f%%)\n", waits / pe_prof.size(), waits * 100.0 / runs);

        auto csl1 = cache_stats(cache_prof_l1);
        auto csl2 = cache_stats(cache_prof_l2);

        printf("\n");
        printf("Cache statistics:\n");
        printf("  L1 (1/1PE, 2048 Units, Line: 64B, Latency: 16):\n");
        printf("    read hit rate:  %11.2f%%\n", get<0>(csl1));
        printf("    read miss:      %12ld (%.3fGB)\n", get<1>(csl1), get<1>(csl1) * 64.0 / 1e9);
        printf("    write hit rate: %11.2f%%\n", get<2>(csl1));
        printf("    write miss:     %12ld (%.3fGB)\n", get<3>(csl1), get<3>(csl1) * 64.0 / 1e9);

        printf("  L2: (1/16PE, 128 Units, Line: 256B, Latency: 28)\n");
        printf("    read hit rate:  %11.2f%%\n", get<0>(csl2));
        printf("    read miss:      %12ld (%.3fGB)\n", get<1>(csl2), get<1>(csl2) * 256.0 / 1e9);
        printf("    write hit rate: %11.2f%%\n", get<2>(csl2));
        printf("    write miss:     %12ld (%.3fGB)\n", get<3>(csl2), get<3>(csl2) * 256.0 / 1e9);
    }

    template <size_t N>
    std::tuple<double, uint64_t, double, uint64_t>
    cache_stats(const std::array<pzcl_profile_cache, N> &profs) const
    {
        size_t read_requests = 0, read_hits = 0;
        size_t write_requests = 0, write_hits = 0;

        for (auto &p : profs) {
            read_requests += p.read_request;
            read_hits += p.read_hit;
            write_requests += p.write_request;
            write_hits += p.write_hit;
        }

        return std::make_tuple(
            read_hits * 100.0 / read_requests,
            read_requests - read_hits,
            write_hits * 100.0 / write_requests,
            write_requests - write_hits);
    }

    pzcl_profile_pe_stats pe_stats;
    pzcl_profile_cache_stats cache_stats_l1;
    pzcl_profile_cache_stats cache_stats_l2;

    std::array<pzcl_profile_pe, 2048> pe_prof;
    std::array<pzcl_profile_cache, 512> cache_prof_l1;
    std::array<pzcl_profile_cache, 128> cache_prof_l2;
};
