import pytest
import tempfile
import os
from scripts.insert_virus import Events

# scripts/test_insert_virus.py


@pytest.fixture
def host_fasta(tmp_path):
    content = ">chr1\n" + "A" * 1000 + "\n"
    file = tmp_path / "host.fa"
    file.write_text(content)
    return str(file)

@pytest.fixture
def virus_fasta(tmp_path):
    content = ">virus1\n" + "T" * 100 + "\n"
    file = tmp_path / "virus.fa"
    file.write_text(content)
    return str(file)

@pytest.fixture
def prob_dict():
    return {
        'p_whole': 0.5,
        'p_rearrange': 0.1,
        'p_delete': 0.1,
        'lambda_split': 2.0,
        'p_overlap': 0.1,
        'p_gap': 0.1,
        'lambda_junction': 2.0,
        'p_host_del': 0.1,
        'lambda_host_del': 5.0
    }

def test_events_init_valid(host_fasta, virus_fasta):
    e = Events(host_fasta, virus_fasta)
    assert hasattr(e, "host")
    assert hasattr(e, "virus")
    assert isinstance(e.host, dict)
    assert isinstance(e.virus, dict)

@pytest.mark.parametrize("seed", [0, -1, "abc", [], None])
def test_events_init_invalid_seed(host_fasta, virus_fasta, seed):
    with pytest.raises(AssertionError):
        Events(host_fasta, virus_fasta, seed=seed)

@pytest.mark.parametrize("min_len,max_len,err", [
    (-1, None, AssertionError),
    (0, None, AssertionError),
    (5, 2, ValueError),
    (None, 1, AssertionError),
])
def test_events_init_invalid_lengths(host_fasta, virus_fasta, min_len, max_len, err):
    with pytest.raises(err):
        Events(host_fasta, virus_fasta, min_len=min_len, max_len=max_len)

def test_events_init_missing_file(tmp_path, virus_fasta):
    with pytest.raises(OSError):
        Events(str(tmp_path / "nofile.fa"), virus_fasta)
    with pytest.raises(OSError):
        Events(virus_fasta, str(tmp_path / "nofile.fa"))


def test__check_probs_missing_key(host_fasta, virus_fasta):
    e = Events(host_fasta, virus_fasta)
    bad_probs = {'p_whole': 0.5}
    with pytest.raises(ValueError):
        e._check_probs(bad_probs)

def test_add_integrations_and_double_add(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_integrations(prob_dict, 2)
    assert hasattr(e, "ints")
    assert len(e.ints) == 2
    with pytest.raises(ValueError):
        e.add_integrations(prob_dict, 1)

def test_add_integrations_too_many(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    # Too many integrations for host size
    with pytest.raises(ValueError):
        e.add_integrations(prob_dict, 10000)

def test_add_integrations_bad_args(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    with pytest.raises(AssertionError):
        e.add_integrations(prob_dict, -1)
    with pytest.raises(AssertionError):
        e.add_integrations(prob_dict, 1, max_attempts=0)
    with pytest.raises(AssertionError):
        e.add_integrations(prob_dict, 1, min_sep=0)

def test_add_episomes_and_double_add(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_episomes(prob_dict, 2)
    assert hasattr(e, "epis")
    assert len(e.epis) == 2
    with pytest.raises(ValueError):
        e.add_episomes(prob_dict, 1)

def test_add_episomes_bad_args(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    with pytest.raises(AssertionError):
        e.add_episomes(prob_dict, -1)
    with pytest.raises(AssertionError):
        e.add_episomes(prob_dict, 1, max_attempts=0)

def test_check_junction_length_raises(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta, min_len=2)
    bad_probs = prob_dict.copy()
    bad_probs['lambda_junction'] = 100
    bad_probs['p_gap'] = 0.5
    bad_probs['p_overlap'] = 0.5
    with pytest.raises(ValueError):
        e.check_junction_length(bad_probs)

def test_save_fasta_and_info(tmp_path, host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_integrations(prob_dict, 1)
    e.add_episomes(prob_dict, 1)
    fasta_file = tmp_path / "out.fa"
    int_info = tmp_path / "int.tsv"
    epi_info = tmp_path / "epi.tsv"
    e.save_fasta(str(fasta_file))
    assert os.path.exists(fasta_file)
    e.save_integrations_info(str(int_info))
    assert os.path.exists(int_info)
    e.save_episomes_info(str(epi_info))
    assert os.path.exists(epi_info)
                                       
def test_add_integrations_creates_ints(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_integrations(prob_dict, 2)
    assert hasattr(e, "ints")

def test_add_integrations_len(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_integrations(prob_dict, 2)
    assert len(e.ints) == 2

def test_add_integrations_double_add_raises(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_integrations(prob_dict, 2)
    with pytest.raises(ValueError):
        e.add_integrations(prob_dict, 1)

def test_add_episomes_creates_epis(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_episomes(prob_dict, 2)
    assert hasattr(e, "epis")

def test_add_episomes_len(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_episomes(prob_dict, 2)
    assert len(e.epis) == 2

def test_add_episomes_double_add_raises(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    e.add_episomes(prob_dict, 2)
    with pytest.raises(ValueError):
        e.add_episomes(prob_dict, 1)

def test_add_episomes_bad_args_negative(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    with pytest.raises(AssertionError):
        e.add_episomes(prob_dict, -1)

def test_add_episomes_bad_args_zero_attempts(host_fasta, virus_fasta, prob_dict):
    e = Events(host_fasta, virus_fasta)
    with pytest.raises(AssertionError):
        e.add_episomes(prob_dict, 1, max_attempts=0)

def test_checkFastaExists_true_host(host_fasta, virus_fasta):
    e = Events(host_fasta, virus_fasta)
    assert e.checkFastaExists(host_fasta)

def test_checkFastaExists_true_virus(host_fasta, virus_fasta):
    e = Events(host_fasta, virus_fasta)
    assert e.checkFastaExists(virus_fasta)

def test_checkFastaExists_false_nonexistent(host_fasta, virus_fasta):
    e = Events(host_fasta, virus_fasta)
    assert not e.checkFastaExists("not_a_file.fa")