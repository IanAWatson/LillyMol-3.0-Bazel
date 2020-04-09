#include <string>

#include "aromatic.h"
#include "substructure.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

//using google::protobuf::TextFormat::ParseFromString;

using testing::UnorderedElementsAre;

class TestSubstructure : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

  protected:
    void _WriteQuery(const char * fname) {
      IWString tmp(fname);
      _query.write_msi(tmp);
    }
};

void
TestSubstructure::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestSubstructure, SingleAtom)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  ASSERT_TRUE(_m.build_from_smiles("C"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, SingleAtomMultipleMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  for (int i = 0; i < 20; ++i) {
    _smiles << "C";
    ASSERT_TRUE(_m.build_from_smiles(_smiles));
    EXPECT_EQ(_query.substructure_search(_m, _sresults), i + 1);
  }
}

TEST_F(TestSubstructure, SingleAtomNot)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        match_as_match: false
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestComposite)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeAnd)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_AND
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeLowPriorityAnd)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_LP_AND
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  EXPECT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
  _smiles = "NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompositeXor)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
    logexp: SS_XOR
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestCompositeEachComponentMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "C"
        }
      }
    }
    match_each_component: 1
    query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "N"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestRespectInitialAtomNumbering)
{
  _string_proto = R"(query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 0
        initial_atom_number: 2
        atom_properties {
          atomic_symbol : "Na"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[Na]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 3);

  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(0));
  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(1));
  EXPECT_EQ(0, e->item(2));
}

TEST_F(TestSubstructure, TestMatchAsMatchNoMatch)
{
  _string_proto = R"(query {
      query_atom {
        match_as_match: false
        id: 0
        atom_properties {
          atomic_symbol : "I"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[I]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestMatchAsMatchMatch)
{
  _string_proto = R"(query {
      query_atom {
        match_as_match: false
        id: 0
        atom_properties {
          atomic_symbol : "I"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[I]C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(e->item(0), 1);

  _smiles = "C[I]";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructure, TestSymmetricMatches)
{
  _string_proto = R"(query {
      perceive_symmetric_equivalents: false
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol : "C"
        }
      }
      query_atom {
        id: 3
        single_bond: 0
        atom_properties {
          atomic_number: 6
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  _query.set_perceive_symmetry_equivalent_matches(true);
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestUniqueEmbeddingsOnly)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "C(F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestIncludeAtomInEmbeddingSmarts)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "[/IWxC](F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(INVALID_ATOM_NUMBER, e->item(0));
  EXPECT_EQ(_m.atomic_number(e->item(1)), 9);
  _query.set_respect_initial_atom_numbering(0);
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 4);
}

TEST_F(TestSubstructure, TestAllVariantsMatch)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      smarts: "*(F)(F)(F)F"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 24);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 5);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestOneEmebeddingPerStartAtom)
{
  // When build from smarts, respect_initial_atom_numbering is in effect.
  _string_proto = R"(query {
      unique_embeddings_only: true
      smarts: "c1ccccc1"
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 6);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestIncludeAtomInEmbeddingProto)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      query_atom {
        id: 0
        include_in_embedding: false
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 9
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 4);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 9);
}

TEST_F(TestSubstructure, TestAtomOrId)
{
  _string_proto = R"(query {
      unique_embeddings_only: true
      query_atom {
        id: 0
        or_id: 75
        atom_properties {
          atomic_symbol : "F"
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FC(F)(F)F";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  for (int i = 0; i < 4; ++i) {
    const Set_of_Atoms * e = _sresults.embedding(0);
    EXPECT_EQ(e->number_elements(), 1);
  }
}

TEST_F(TestSubstructure, TestInitialAtomNumberSparse)
{
  _string_proto = R"(query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 0
        initial_atom_number: 8
        atom_properties {
          atomic_symbol : "O"
        }
      }
      query_atom {
        id: 1
        initial_atom_number: 2
        atom_properties {
          atomic_symbol : "C"
        }
        double_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC(=O)C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 9);
  EXPECT_EQ(e->item(0), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(1), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(2), 1);
  EXPECT_EQ(e->item(3), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(4), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(5), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(6), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(7), INVALID_ATOM_NUMBER);
  EXPECT_EQ(e->item(8), 2);
}

TEST_F(TestSubstructure, TestRingIdMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        ring_id: 8
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        ring_id: 8
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
        aromatic_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 2);
}

TEST_F(TestSubstructure, TestRingIdNotMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        ring_id: 2
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        ring_id: 8
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
        aromatic_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestFusedSystemIdMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fused_system_id: 3
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        fused_system_id: 3
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_symbol : "C"
          aromatic: false
        }
        single_bond: 1
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1nnc2cn(C)ccc12";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->number_elements(), 3);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructure, TestFusedSystemIdNoMatch)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fused_system_id: 3
        atom_properties {
          atomic_symbol: "C"
          aromatic: true
        }
      }
      query_atom {
        id: 1
        fused_system_id: 9
        atom_properties {
          atomic_symbol : "N"
          aromatic: true
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1ccccc1c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 6);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
  EXPECT_EQ(_m.atomic_number(e->item(1)), 7);
}

TEST_F(TestSubstructure, TestFragmentIdMatchesSameFrag)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 3
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestFragmentIdNoMatches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 9
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestFragmentIdMatchesDifferentFrags)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        fragment_id: 3
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        fragment_id: 9
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C.N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsSingle)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsMultipleAttach)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0"
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
          aromatic: false
        }
        bond_smarts: "- 0 1"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestBondSmartsDouble)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsTriple)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "# 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C#NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsArom)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: ": 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1nccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestBondSmartsAny)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "~ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C-N.C=N.C#N.c1ncccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 5);
}

TEST_F(TestSubstructure, TestBondSmartsRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C-N.C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}


TEST_F(TestSubstructure, TestBondSmartsDoubleAndRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "@= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsDoubleNotRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "!@= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N=C.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestBondSmartsSingleOrDouble)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "-,= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C=N-C.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestBondSmartsNonRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "!@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();


  _smiles = "C-N.C1NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestBondSmartsXor)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "@^= 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestBondSmartsSingleDoubleAndRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        bond_smarts: "-,=;@ 0"
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
}

TEST_F(TestSubstructure, TestQueryBondSimple)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        query_bond {
          bond_type: SS_SINGLE_BOND
          other_end: 0
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C1NC1.C1C=NC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestQueryBondComplex)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        query_bond {
          bond_type: SS_SINGLE_BOND
          bond_type: SS_TRIPLE_BOND
          other_end: 0
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC#N";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
}

TEST_F(TestSubstructure, TestPreferenceDefault)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 5
        }
        preference {
          ncon: 1
          preference_value: 1
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  // By default, low preference hit are removed.

  _smiles = "NC.NCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);

  _smiles = "NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestPreferenceSortPreferences)
{
  _string_proto = R"(query {
      sort_by_preference_value: true
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 5
        }
        preference {
          ncon: 1
          preference_value: 1
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  // By default, low preference hit are removed.

  _smiles = "NC.NCC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);

  _smiles = "NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestPreferenceSumPreferences)
{
  _string_proto = R"(query {
      sort_by_preference_value: true
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
        preference {
          ncon: 2
          preference_value: 1
        }
        preference {
          ncon: 1
          preference_value: 5
        }
        preference {
          ring_bond_count: 2
          preference_value: 6
        }
        sum_all_preference_hits: true
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "NC.NCC.N1CC1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
  EXPECT_EQ(_m.ring_bond_count(e->item(0)), 2);

  _smiles = "N1CC1.NCC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);
  e = _sresults.embedding(0);
  EXPECT_EQ(_m.ncon(e->item(0)), 2);
  EXPECT_EQ(_m.ring_bond_count(e->item(0)), 2);
}

TEST_F(TestSubstructure, TestElementsNeededMatches)
{
  _string_proto = R"(query {
      elements_needed {
        atomic_number: 9
        hits_needed: 1
      }
      elements_needed {
        atomic_number: 6
        min_hits_needed: 1
      }
      elements_needed {
        atomic_number: 7
        max_hits_needed: 1
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestElementsNeededNoMatches)
{
  _string_proto = R"(query {
      elements_needed {
        atomic_number: 9
        hits_needed: 2
      }
      elements_needed {
        atomic_number: 6
        min_hits_needed: 1
      }
      elements_needed {
        atomic_number: 7
        max_hits_needed: 1
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestElementHitsNeededMatches)
{
  _string_proto = R"(query {
      element_hits_needed {
        atomic_number: 9
        hits_needed: 0
      }
      elements_needed {
        atomic_number: 6
        min_hits_needed: 1
      }
      elements_needed {
        atomic_number: 7
        max_hits_needed: 1
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestElementHitsNeededNoMatches)
{
  _string_proto = R"(query {
      element_hits_needed {
        atomic_number: 9
        hits_needed: 1
      }
      elements_needed {
        atomic_number: 6
        min_hits_needed: 1
      }
      elements_needed {
        atomic_number: 7
        max_hits_needed: 1
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 7
        }
        single_bond: 0
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "FCN";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestChiral4Matches)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 7
        }
        single_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 1
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 1
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 6
        }
        single_bond: 4
      }
      chiral_centre {
        center: 1
        top_front {
          atom_number: 0
        }
        top_back {
          atom_number: 2
        }
        left_down {
          atom_number: 3
        }
        right_down {
          atom_number: 4
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C[C@@](F)(N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "C[C@](F)(N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructure, TestChiral3Matches)
{
  _string_proto = R"(query {
      smarts: "C[CH](N)CC";
      chiral_centre {
        center: 1
        top_front {
          atom_number: 0
        }
        top_back {
          atom_number: 2
        }
        left_down {
          h_or_lp: "H"
        }
        right_down {
          atom_number: 3
        }
      }
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "C[C@@H](N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  _smiles = "C[C@H](N)CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);
}


TEST_F(TestSubstructure, TestAnyLength)
{
  _string_proto = R"(query {
      smarts: "C#[#{hello}D2]=[#{world}D>1]";
    }
  )";

  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "[Cr][world]=[hello]#CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructure, TestCompsite)
{
  _string_proto = R"(query {
      smarts: "CC";
    }
    query {
      smarts: "NN";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1cccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(2, 3));
  }

  _smiles = "CC.c1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(0, 1));
  }
}

TEST_F(TestSubstructure, TestCompsiteXor)
{
  _string_proto = R"(query {
      smarts: "CC";
    }
    logexp: SS_XOR
    query {
      smarts: "NN";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "c1cccnc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.CC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "NN.C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_EQ(_sresults.number_embeddings(), 2);

  for (const auto* e : _sresults.embeddings())
  {
    EXPECT_EQ(e->number_elements(), 2);
    EXPECT_THAT(*e, UnorderedElementsAre(0, 1));
  }

  _smiles = "CC.NC";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
  // A surprising result, no embeddings. This is correct because
  // the last query done, looking for NN, failed. But at that point
  // the result of the XOR was available. But at that time,
  // _sresults was holding a non match.
  EXPECT_EQ(_sresults.number_embeddings(), 0);
}

TEST_F(TestSubstructure, TestManyMatches)
{
  _string_proto = R"(query {
      smarts: "CC(C)(C)c1ccc(C(C)(C)C)cc1";
    }
  )";

  SubstructureSearch::SubstructureQuery proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &proto));

  ASSERT_TRUE(_query.ConstructFromProto(proto)) << "Cannot parse proto " << proto.ShortDebugString();

  _smiles = "CC(C)(C)c1ccc(C(C)(C)C)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 144);
}

}  // namespace
