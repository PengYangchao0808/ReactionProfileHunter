from pathlib import Path


def test_gaussian_route_block_wraps_without_splitting_tokens():
    from rph_core.utils.qc_interface import _format_gaussian_route_block

    route = (
        "B3LYP/def2SVP em=GD3BJ SCRF=(PCM,Solvent=Acetone) Opt=CalcFC Freq "
        "Pop=NBO"
    )
    block = _format_gaussian_route_block(route, max_cols=80)
    lines = block.splitlines()

    assert lines
    assert lines[0].startswith("#p ")
    for ln in lines:
        assert len(ln) <= 80

    # Regression: never split 'Pop' into 'P' + 'op=...'
    for i in range(len(lines) - 1):
        assert not (lines[i].endswith(" P") and lines[i + 1].lstrip().startswith("op="))


def test_gaussian_interface_write_input_file_wraps_route(tmp_path: Path):
    from rph_core.utils.qc_interface import GaussianInterface

    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    out = tmp_path / "out"
    out.mkdir()
    gjf = out / "m.gjf"

    gi = GaussianInterface(config={"executables": {"gaussian": {"use_wrapper": True}}})
    gi.write_input_file(
        xyz_file=xyz,
        gjf_file=gjf,
        route=(
            "#p B3LYP/def2SVP em=GD3BJ SCRF=(PCM,Solvent=Acetone) "
            "Opt=(CalcFC,NoEigenTest) Freq Pop=NBO"
        ),
        title="t",
    )

    text = gjf.read_text().splitlines()
    route_start = next(i for i, ln in enumerate(text) if ln.startswith("#p "))

    route_lines = []
    for ln in text[route_start:]:
        if not ln.strip():
            break
        route_lines.append(ln)

    assert route_lines
    for ln in route_lines:
        assert len(ln) <= 80
    assert any(ln.startswith(" ") for ln in route_lines[1:])


def test_gaussian_interface_injects_nbo_keylist_when_requested(tmp_path: Path):
    from rph_core.utils.qc_interface import GaussianInterface

    xyz = tmp_path / "m.xyz"
    xyz.write_text("1\ncomment\nH 0.0 0.0 0.0\n")

    out = tmp_path / "out"
    out.mkdir()
    gjf = out / "m.gjf"

    gi = GaussianInterface(
        config={
            "executables": {"gaussian": {"use_wrapper": True}},
            "theory": {"optimization": {"nbo_keylist": "archive\n"}},
        }
    )
    gi.write_input_file(
        xyz_file=xyz,
        gjf_file=gjf,
        route="#p B3LYP/def2SVP Opt Freq Pop=(NBORead)",
        title="t",
    )

    text = gjf.read_text()
    assert "$NBO" in text
    assert "archive" in text
    assert "$END" in text
