cc_binary(
    name = "concat_files",
    srcs = [
        "concat_files.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "correlate",
    srcs = [
       "correlate.cc",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:iw_tabular_data",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "fetch_smiles_quick",
    srcs = [
        "fetch_smiles_quick.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwcut",
    srcs = [
        "iwcut.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "unique_rows",
    srcs = [
        "unique_rows.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)
