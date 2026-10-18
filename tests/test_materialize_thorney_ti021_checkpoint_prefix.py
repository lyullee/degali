from tools.materialize_thorney_ti021_checkpoint_prefix import PREFIX_SCHEMA


def test_checkpoint_prefix_schema_is_versioned():
    assert PREFIX_SCHEMA == "degali.thorney-ti021-checkpoint-prefix.v1"
