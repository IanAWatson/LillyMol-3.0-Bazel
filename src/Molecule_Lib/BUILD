package(default_visibility = ["//visibility:public"])

#genrule(
#  name = 'gfortran_iw',
#  srcs = [
#    'u3b.f',
#  ],
#  outs = ['u3b.o'],
#  cmd = 'gfortran -O3 -c -o $(OUTS) u3b.f'
#)

proto_library(
    name = "substructure_proto",
    srcs = [
        "substructure.proto"
    ],
)

cc_proto_library(
    name = "substructure_cc_proto",
    deps = [
        ":substructure_proto",
    ],
)

cc_library(
    name = "iwmolecule",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "allowed_elements.cc",
        "aromatic.cc",
        "atom.cc",
        "atom_alias.cc",
        "atom_typing.cc",
        "bond.cc",
        "bond_list.cc",
        "careful_frag.cc",
        "charge_assigner.cc",
        "charge_calculation.cc",
        "charmm.cc",
        "chiral_centre.cc",
        "cif.cc",
        "cis_trans_bond.cc",
        "coordinates.cc",
        "csubstructure.cc",
        "dihedral.cc",
        "donor_acceptor.cc",
        "element.cc",
        "element_hits_needed.cc",
        "ematch.cc",
        "etrans.cc",
        "frag.cc",
        "ISIS_Atom_List.cc",
        "is_actually_chiral.cc",
        "iwrcb.cc",
        "iwqry_wstats.cc",
        "iwrnm.cc",
        "iwsubstructure.cc",
        "mdl_atom_record.cc",
        "mdl_file_data.cc",
        "mdl_molecule.cc",
        "mdl_v30.cc",
        "misc2.cc",
        "molecule.cc",
        "molecule_to_query.cc",
        "moleculeb.cc",
        "moleculed.cc",
        "moleculeh.cc",
        "moleculer.cc",
        "molecule_smarts.cc",
        "mrk.cc",
        "msi.cc",
        "numass.cc",
        "ostream_and_type.cc",
        "output.cc",
        "parse_smarts_tmp.cc",
        "path.cc",
        "path_scoring.cc",
        "pdb.cc",
        "pearlman.cc",
        "reaction_match_conditions.cc",
        "ring_bond_iterator.cc",
        "rmele.cc",
        "rotbond_common.cc",
        "rwsubstructure.cc",
        "rwsubstructure_proto.cc",
        "set_of_atoms.cc",
        "smi.cc",
        "smiles.cc",
        "smiles_support.cc",
        "ss_atom_env.cc",
        "ss_bonds.cc",
        "ss_ring.cc",
        "ss_ring_base.cc",
        "ss_ring_sys.cc",
        "substructure_chiral.cc",
        "substructure_results.cc",
        "symm_class_can_rank.cc",
        "standardise.cc",
        "substructure_a.cc",
        "substructure_nmab.cc",
        "substructure_spec.cc",
        "substructure_env.cc",
        "symmetry.cc",
        "target.cc",
        "temp_detach_atoms.cc",
        "toggle_kekule_form.cc",
        "tokenise_atomic_smarts.cc",
        "unique.cc",
        "u3b.c",
        "_substructure.cc",
        "_qry_wstats.cc",
    ],
    hdrs = [
        "allowed_elements.h",
        "aromatic.h",
        "atom.h",
        "atom_alias.h",
        "atom_typing.h",
        "bond.h",
        "bond_list.h",
        "careful_frag.h",
        "charge_assigner.h",
        "charge_calculation.h",
        "cis_trans_bond.h",
        "coordinates.h",
        "chiral_centre.h",
        "collection_template.h",
        "dihedral.h",
        "donor_acceptor.h",
        "element.h",
        "ematch.h",
        "etrans.h",
        "istream_and_type.h",
        "is_actually_chiral.h",
        "iwmtypes.h",
        "iwrcb.h",
        "iwrnm.h",
        "iwstandard.h",
        "marvin.h",
        "mdl.h",
        "mdl_atom_record.h",
        "mdl_file_data.h",
        "mdl_molecule.h",
        "misc2.h",
        "molecule.h",
        "molecule_to_query.h",
        "molecule_tripos.h",
        "moleculed.h",
        "moleculeh.h",
        "moleculej.h",
        "molecule_arom.h",
        "molecule_cif.h",
        "molecule_ctb.h",
        "molecule_main.h",
        "molecule_marvin.h",
        "molecule_mdl.h",
        "molecule_smi.h",
        "molecule_smarts.h",
        "numass.h",
        "ostream_and_type.h",
        "output.h",
        "path.h",
        "path_scoring.h",
        "parse_smarts_tmp.h",
        "pearlman.h",
        "qry_wstats.h",
        "readmdl.h",
        "reaction_match_conditions.h",
        "rotbond_common.h",
        "rmele.h",
        "rwmolecule.h",
        "rwsubstructure.h",
        "set_of_atoms.h",
        "smiles.h",
        "space_vector.h",
        "substructure.h",
        "symmetry.h",
        "target.h",
        "temp_detach_atoms.h",
        "tmpsssr.h",
        "toggle_kekule_form.h",
        "tokenise_atomic_smarts.h",
    ],
    deps = [
        ":substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iw_tdt:iw_tdt",
        "//Foundational/mtrand:iwrandom",
        "//Foundational/xmlParser:xmlparser",
    ],
)

cc_library(
    name = "iwsubstructure",
    srcs = [
    ],
    hdrs = [
    ],
    deps = [
    ],
)

cc_library(
    name = "moleculeio",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "charmm.cc",
        "cif.cc",
        "inchi_dummy.cc",
        "marvin.cc",
        "mdl.cc",
        "mdl_v30.cc",
        "mmod.cc",
        "rwmolecule.cc",
        "tripos.cc",
    ],
    hdrs = [
        "tripos.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)

cc_library(
    name = "iwreaction",
    srcs = [
       "iwreaction.cc",
       "iwreaction2.cc",
        "rxn_3d_replace.cc",
        "rxn_bond_angle.cc",
        "rxn_bond_length.cc",
        "rxn_dihedral_angle.cc",
        "rxn_matched_sors.cc",
    ],
    hdrs = [
        "iwreaction.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)

cc_library(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint.cc",
    ],
    hdrs = [
        "atom_pair_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "ec_fingerprint",
    srcs = [
        "ec_fingerprint.cc",
    ],
    hdrs = [
        "ec_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint.cc",
    ],
    hdrs = [
        "linear_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_test(
    name = "linear_fingerprint_test",
    srcs = ["linear_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":linear_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "ec_fingerprint_test",
    srcs = ["ec_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":ec_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "atom_pair_fingerprint_test",
    srcs = ["atom_pair_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":atom_pair_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_test",
    srcs = ["substructure_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_spec_test",
    srcs = ["substructure_spec_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_single_query_test",
    srcs = ["substructure_single_query_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_env_test",
    srcs = ["substructure_env_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_ring_test",
    srcs = ["substructure_ring_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_ring_system_test",
    srcs = ["substructure_ring_system_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_smarts_test",
    srcs = ["substructure_smarts_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "standardise_test",
    srcs = ["standardise_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)
