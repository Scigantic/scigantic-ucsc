import os

import pytest

from scigantic_ucsc import toga

ROW = ("ENST00000308824.PCMTD2.47\t1000\t-\t183903\t197236\t0,0,200\t5\t377,124,\t0,2655,"
       "\tENST00000308824.PCMTD2\tchr20:64255747-64276226\tCZRN01000001:183903-197236\t0.99\t343\t0.30"
       "\t0.03\t0.09\t0.99\t0.63\tIntact\t0.99\t0.98\t0.97\t0.0\t1\t1\t<TT>big html blob")


def test_parse_entry():
    g = toga.parse_entry("CZRN01000001", 183903, 197236, ROW)
    assert (g.status, g.strand, g.gene_symbol) == ("Intact", "-", "PCMTD2")
    assert g.fields["chain_score"] == 0.99 and g.fields["perc_intact_ign_missing"] == 0.99
    assert "prot_alignment" not in g.fields and "<TT>" not in str(g.fields)


def test_accession_paths():
    assert toga._prefix("GCF_001458135.1") == "hubs/GCF/001/458/135/GCF_001458135.1/bbi/"
    assert toga._prefix("GCA_009769595.1").startswith("hubs/GCA/009/769/595/")
    with pytest.raises(ValueError):
        toga._prefix("marmoset")


@pytest.mark.skipif(not os.environ.get("SCIGANTIC_UCSC_LIVE"), reason="live bucket")
def test_live_marmoset():
    genes = toga.toga_genes("GCF_001458135.1", "CZRN01000001", 0, 5_000_000)
    assert len(genes) > 100
    assert all(g.status for g in genes)
    intact = toga.toga_genes("GCF_001458135.1", "CZRN01000001", 0, 5_000_000, status="intact")
    assert 0 < len(intact) <= len(genes)
    with pytest.raises(ValueError, match="no TOGA"):
        toga.toga_key("GCF_000000000.1")


def test_open_remote_times_out_even_when_the_open_blocks_the_gil(monkeypatch):
    import time

    from scigantic_ucsc import _http

    # A busy loop holds the GIL the way pyBigWig.open does; a thread-based
    # timeout would never fire here.
    monkeypatch.setattr(_http, "_PROBE", "while True: pass\n")
    t = time.time()
    with pytest.raises(TimeoutError, match="took longer"):
        _http.open_remote("never.bb", timeout=1.0)
    assert time.time() - t < 10
