cc_binary(
    name = "gfp_leader",
    srcs = [
        "gfp_leader_v2.cc",
        "leader.h",
    ],
    deps = [
        "gfp",
        ":sparse_collection",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_spread",
    srcs = [
        "gfp_spread_v2.cc",
    ],
    deps = [
        ":gfp",
        ":spread_v2",
        ":smiles_id_dist",
        ":sparse_collection",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_lnearneighbours",
    srcs = [
        "gfp_lnearneighbours.cc",
        "fp_and_smiles.cc",
    ],
    deps = [
        ":gfp",
        ":neighbour_list",
        ":sparse_collection",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "nplotnn",
    srcs = [
        "nplotnn.cc",
    ],
    deps = [
        ":distance_scaling",
        ":smiles_id_dist",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_distance_matrix",
    srcs = [
        "gfp_distance_matrix.cc",
    ],
    deps = [
        ":gfp",
        ":build_pool",
        ":sparse_collection",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_library(
    name = "spread_v2",
    srcs = [
        "spread_v2.cc",
    ],
    hdrs = [
        "spread_v2.h",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "build_pool",
    srcs = [
        "build_pool.cc",
    ],
    hdrs = [
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "smiles_id_dist",
    srcs = [
        "smiles_id_dist.cc",
    ],
    hdrs = [
        "smiles_id_dist.h",
    ],
    deps = [
        ":gfp",
    ],
)


cc_library(
    name = "gfp",
    srcs = [
        "bic_table.cc",
        "gfp.cc",
        "dyfp.cc",
        "fixed_size_counted_fingerprint_uchar.cc",
        "_fixed_size_counted_fingerprint_uchar.cc",
        "_fixed_size_counted_fingerprint_uint.cc",
        "multi_conformer.cc",
        "sparsefp.cc",
        "tversky.cc",
        "various_distance_metrics.cc",
    ],
    hdrs = [
        "gfp.h",
        "dyfp.h",
        "fixed_size_counted_fingerprint.h",
        "multi_conformer.h",
        "sparsefp.h",
        "sparse_collection.h",
        "tversky.h",
        "various_distance_metrics.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iw_tdt:iw_tdt",
    ],
)

cc_library(
    name = "neighbour_list",
    srcs = [
        "neighbour_list.cc",
        "_neighbour_list_fp_and_smiles.cc",
    ],
    hdrs = [
        "neighbour_list.h",
    ],
    deps = [
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "sparse_collection",
    srcs = [
        "sparse_collection.cc",
    ],
    hdrs = [
        "sparse_collection.h",
    ],
    deps = [
        ":gfp",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_library(
    name = "distance_scaling",
    srcs = [
        "distance_scaling.cc",
    ],
    hdrs = [
        "distance_scaling.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwmisc",
    ],
)
