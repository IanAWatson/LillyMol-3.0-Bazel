cc_binary(
    name = "fileconv",
    srcs = [
        "fileconv.cc",
        "fix_structures.cc",
        "fix_structures.h",
        "do_remove_duplicate_fragments.cc",
        "do_remove_duplicate_fragments.h",
        "known_fragment_data.cc",
        "known_fragment_data.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tsubstructure",
    srcs = [
        "tsubstructure.cc",
        "tsubstructure_fp.cc",
        "tsubstructure_fp.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "unique_molecules",
    srcs = [
        "unique_molecules.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "msort",
    srcs = [
        "msort.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "mol2qry",
    srcs = [
        "mol2qry.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rotatable_bonds",
    srcs = [
        "rotatable_bonds.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tp_first_pass",
    srcs = [
        "tp_first_pass.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz"
    ],
)

cc_binary(
    name = "iwdemerit",
    srcs = [
        "iwdemerit.cc",
        "substructure_demerits.cc",
        "substructure_demerits.h",
        "demerit.cc",
        "demerit.h",
        "qry_and_demerit.h",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz"
    ],
)

cc_binary(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:atom_pair_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ec_fingerprint",
    srcs = [
        "ec_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:ec_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:linear_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)
