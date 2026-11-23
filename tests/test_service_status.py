import json
import subprocess
import sys
from unittest.mock import Mock

import pytest

from pisa_utils import run
from pisa_utils.models.validation_error_handlers import ErrorLabels
from pisa_utils.run_pisa import StatusReporter


PARSERS = (
    "ConvertInterfaceXMLToJSONs",
    "ConvertAssemblyXMLToJSON",
    "CompileInterfaceSummaryJSON",
    "ConvertAssemblyListToJSON",
    "ConvertInterfaceListToJSON",
    "ConvertComponentsListToJSON",
    "PostProcessComplexTable",
    "PostProcessInterfaceDetailsList",
    "PostProcessInterfaceTable",
    "PostProcessInterfaceComponentResidueTable",
    "PostProcessInterfacingBondsTable",
)


def test_status_report_uses_custom_filename(tmp_path):
    report = StatusReporter(str(tmp_path), filename="result.json")
    report.write_service_status(0, "Success")

    assert json.loads((tmp_path / "result.json").read_text()) == {
        "status_code": 0,
        "message": "Success",
    }
    assert not (tmp_path / "status.json").exists()


@pytest.fixture
def service_setup(monkeypatch, tmp_path):
    output_dir = tmp_path / "json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pisa_analysis_service",
            "--input_cif",
            "input.cif",
            "--output_xml",
            str(tmp_path / "xml"),
            "--output_json",
            str(output_dir),
            "--compress_output",
        ],
    )
    monkeypatch.setattr(run, "validate_args", lambda args: args)
    pisa = Mock(
        return_value=("assembly", "interfaces", "assemblies", "ifaces", "monomers")
    )
    monkeypatch.setattr(run, "run_pisa_service", pisa)
    parsers = []
    for name in PARSERS:
        parser = Mock()
        parsers.append(parser)
        monkeypatch.setattr(run, name, Mock(return_value=parser))
    return output_dir, pisa, parsers


def test_success_status_is_published_after_all_parsers(service_setup):
    output_dir, pisa, parsers = service_setup
    output_dir.mkdir()
    (output_dir / "status.json").write_text('{"status_code": 1}')

    def check_status_absent():
        assert not (output_dir / "status.json").exists()

    pisa.side_effect = lambda **kwargs: (
        check_status_absent()
        or ("assembly", "interfaces", "assemblies", "ifaces", "monomers")
    )
    for parser in parsers:
        parser.parse.side_effect = check_status_absent

    run.service()

    assert json.loads((output_dir / "status.json").read_text()) == {
        "status_code": 0,
        "message": "Success",
    }
    assert all(parser.parse.call_count == 1 for parser in parsers)
    assert list(output_dir.iterdir()) == [output_dir / "status.json"]


@pytest.mark.parametrize("failure", ["pisa", "parser"])
def test_failure_status_and_exception_are_preserved(service_setup, failure):
    output_dir, pisa, parsers = service_setup
    error = RuntimeError(f"{failure} failed")
    if failure == "pisa":
        pisa.side_effect = error
    else:
        parsers[1].parse.side_effect = error

    with pytest.raises(RuntimeError, match=f"{failure} failed") as raised:
        run.service()

    assert raised.value is error
    assert json.loads((output_dir / "status.json").read_text()) == {
        "status_code": 999,
        "message": f"RuntimeError: {failure} failed",
        "type": ErrorLabels.PARSING.value,
    }
    if failure == "parser":
        parsers[0].parse.assert_called_once()
        parsers[2].parse.assert_not_called()
    else:
        assert all(not parser.parse.called for parser in parsers)


@pytest.mark.parametrize(
    ("returncode", "error_label"),
    [(7, ErrorLabels.PISA_NO_MOL_REF_DATA), (-9, ErrorLabels.PISA_UNKNOWN)],
)
def test_subprocess_failure_preserves_returncode(
    service_setup, returncode, error_label
):
    output_dir, pisa, parsers = service_setup
    error = subprocess.CalledProcessError(returncode, ["pisa", "-analyse"])
    pisa.side_effect = error

    with pytest.raises(subprocess.CalledProcessError) as raised:
        run.service()

    assert raised.value is error
    assert json.loads((output_dir / "status.json").read_text()) == {
        "status_code": returncode,
        "message": f"CalledProcessError: {error}",
        "type": error_label.value,
    }
    assert all(not parser.parse.called for parser in parsers)
