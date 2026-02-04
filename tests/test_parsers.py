import json
from pathlib import Path
from unittest import TestCase

from pisa_utils.parsers import (
    CompileInterfaceSummaryJSON,
    ConvertAssemblyListToJSON,
    ConvertAssemblyXMLToJSON,
    ConvertComponentsListToJSON,
    ConvertInterfaceListToJSON,
    ConvertInterfaceXMLToJSONs,
)

# TODO add tests for PDB files


def remove_files(path: Path):
    """
    For setup() method in tests - remove all files in a directory, or create

    :param path: Path to directory. Directory itself not removed.
    :type path: Path
    """
    if path.exists() and path.is_dir():
        for file in path.iterdir():
            if file.is_file():
                file.unlink()
    else:
        path.mkdir(parents=True, exist_ok=True)


class TestConvertAssemblyXMLToJSON(TestCase):
    """
    Tests for ConvertAssemblyXMLToJSON class.
    """

    def setUp(self):
        super().setUp()

        self.maxDiff = None

        self.base_input_dir = Path("tests/data/")

        # Remove any existing output data
        output_path = Path("tests/data/actual_output/")
        output_path.mkdir(parents=True, exist_ok=True)
        remove_files(output_path)

    def test_parse_multi_assembly_xml(self):

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_assembly_multi_asmset.xml"
        )
        self.output_json = self.base_input_dir.joinpath(
            "actual_output", "3hax_assembly_multi_asmset.json"
        )
        self.expected_json = self.base_input_dir.joinpath(
            "expected_output", "3hax_assembly_multi_asmset.json"
        )

        # Run
        self.converter = ConvertAssemblyXMLToJSON(
            path_xml=self.input_xml,
            path_json=self.output_json,
            path_interface_jsons=self.base_input_dir.joinpath(
                "expected_output",
                "interfaces",
                "3hax_interfaces",
            ),
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        expected = self.expected_json.read_text().strip()
        actual = self.output_json.read_text().strip()

        self.assertEqual(
            expected,
            actual,
            msg="Assembly XML->JSON not parsed correctly for multiple assembly sets.",
        )

    def test_parse_single_assembly_xml(self):

        self.maxDiff = None

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_assembly_single_asmset.xml"
        )
        self.output_json = self.base_input_dir.joinpath(
            "actual_output", "3hax_assembly_single_asmset.json"
        )
        self.expected_json = self.base_input_dir.joinpath(
            "expected_output", "3hax_assembly_single_asmset.json"
        )

        # Run
        self.converter = ConvertAssemblyXMLToJSON(
            path_xml=self.input_xml,
            path_json=self.output_json,
            path_interface_jsons=self.base_input_dir.joinpath(
                "expected_output",
                "interfaces",
                "3hax_interfaces",
            ),
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        expected = self.expected_json.read_text().strip()
        actual = self.output_json.read_text().strip()

        self.assertEqual(
            expected,
            actual,
            msg="Assembly XML->JSON not parsed correctly for a single assembly set.",
        )

    def test_parse_assembly_xml_no_asmset(self):
        """
        Test parsing of assembly XML files with no assemblies defined.
        """

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_assembly_no_asmset.xml"
        )
        self.output_json = self.base_input_dir.joinpath(
            "actual_output", "3hax_assembly_no_asmset.json"
        )
        self.expected_json = self.base_input_dir.joinpath(
            "expected_output", "3hax_assembly_no_asmset.json"
        )

        # Run
        self.converter = ConvertAssemblyXMLToJSON(
            path_xml=self.input_xml,
            path_json=self.output_json,
            path_interface_jsons=self.base_input_dir.joinpath(
                "expected_output",
                "interfaces",
                "3hax_interfaces",
            ),
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        expected = self.expected_json.read_text().strip()
        actual = self.output_json.read_text().strip()

        self.assertEqual(
            expected,
            actual,
            msg="Assembly XML->JSON not parsed correctly for no assemblies defined.",
        )


class TestConvertInterfaceXML(TestCase):
    """
    Tests for ConvertInterfaceXMLToJSON class.
    """

    def setUp(self):
        super().setUp()

        self.maxDiff = None

        self.base_input_dir = Path("tests/data/")

        # Remove any existing output data
        output_path = Path("tests/data/actual_output/")
        output_path.mkdir(parents=True, exist_ok=True)
        remove_files(output_path)

    def test_parse_multi_interface_xml(self):

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_interfaces_multi.xml"
        )
        self.output_json_dir = self.base_input_dir.joinpath(
            "actual_output", "interfaces", "3hax_interfaces_multi"
        )
        self.expected_json_dir = self.base_input_dir.joinpath(
            "expected_output", "interfaces", "3hax_interfaces_multi"
        )

        # Run
        self.converter = ConvertInterfaceXMLToJSONs(
            path_xml=self.input_xml,
            path_jsons=self.output_json_dir,
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        for i in (1, 2, 3):
            expected = (
                self.expected_json_dir.joinpath(f"interface_{i}.json")
                .read_text()
                .strip()
            )

            actual = (
                self.output_json_dir.joinpath(f"interface_{i}.json").read_text().strip()
            )

            self.assertEqual(
                expected,
                actual,
                msg="Interface XML->JSON not parsed correctly for interface " f"{i}.",
            )

    def test_parse_single_interface_xml(self):

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_interfaces_single.xml"
        )
        self.output_json_dir = self.base_input_dir.joinpath(
            "actual_output", "interfaces", "3hax_interfaces_single"
        )
        self.expected_json_dir = self.base_input_dir.joinpath(
            "expected_output", "interfaces", "3hax_interfaces_single"
        )

        # Run
        self.converter = ConvertInterfaceXMLToJSONs(
            path_xml=self.input_xml,
            path_jsons=self.output_json_dir,
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        expected = (
            self.expected_json_dir.joinpath("interface_1.json").read_text().strip()
        )

        actual = self.output_json_dir.joinpath("interface_1.json").read_text().strip()

        self.assertEqual(
            expected,
            actual,
            msg="Interface XML->JSON not parsed correctly for single interface.",
        )

    def test_parse_no_interface_xml(self):
        """
        Test parsing of interface XML files with no interfaces defined.
        """

        self.input_xml = self.base_input_dir.joinpath(
            "mock_data", "3hax_interfaces_none.xml"
        )
        self.output_json_dir = self.base_input_dir.joinpath(
            "actual_output", "interfaces", "3hax_interfaces_none"
        )

        # Run
        self.converter = ConvertInterfaceXMLToJSONs(
            path_xml=self.input_xml,
            path_jsons=self.output_json_dir,
            path_structure_file=str(
                self.base_input_dir.joinpath("mock_data", "3hax.cif")
            ),
        )
        self.converter.parse()

        # Check
        actual_files = list(self.output_json_dir.iterdir())

        self.assertEqual(
            0,
            len(actual_files),
            msg="Interface XML->JSON parsed interfaces when none were defined.",
        )


class TestConvertInterfaceListToJSON(TestCase):
    """
    Tests for ConvertInterfaceListToJSON class.
    """

    def setUp(self):
        super().setUp()

        self.maxDiff = None

        self.base_data_dir = Path("tests/data/")
        self.base_expected_dir = self.base_data_dir.joinpath(
            "expected_output/list_results/"
        )
        self.base_output_dir = self.base_data_dir.joinpath(
            "actual_output/list_results/"
        )
        self.base_input_dir = self.base_data_dir.joinpath("mock_data/list_results/")

        # Remove any existing output data
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        remove_files(self.base_output_dir)

    def test_parse_arbitrary_minimal(self):
        """
        Arbitrary minimal interface list XML file.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_minimal.json")
        converter = ConvertInterfaceListToJSON(
            path_txt=str(
                self.base_input_dir.joinpath(
                    "arbitrary_data_minimal", "interfaces_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_minimal", "interfaces_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary minimal interface list XML->JSON not parsed correctly.",
        )

    def test_parse_6nxr(self):
        """
        Test interface list XML file from PISA analysis of 6nxr.
        """
        path_output = self.base_output_dir.joinpath("interfaces_extended.json")
        # Run
        converter = ConvertInterfaceListToJSON(
            path_txt=str(
                self.base_input_dir.joinpath("6nxr", "interfaces_extended.txt")
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr", "interfaces_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="6nxr interface list XML->JSON not parsed correctly.",
        )

    def test_parse_variable_column_widths(self):
        """
        Test interface list XML file with variable column widths.
        """
        path_output = self.base_output_dir.joinpath(
            "interfaces_extended_variable_widths.json"
        )

        # Run
        converter = ConvertInterfaceListToJSON(
            path_txt=str(
                self.base_input_dir.joinpath(
                    "6nxr_variable_column_widths", "interfaces_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr", "interfaces_extended.json"
        )

        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Variable column widths interface list XML->JSON not parsed correctly.",
        )

    def test_parse_arbitrary_maximal(self):
        """
        Arbitrary maximal complexity interface list XML file. Should cover all edge
        cases.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_maximal.json")
        converter = ConvertInterfaceListToJSON(
            path_txt=str(
                self.base_input_dir.joinpath(
                    "arbitrary_data_maximal", "interfaces_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_maximal", "interfaces_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary maximal interface list XML->JSON not parsed correctly.",
        )

    def test_parse_no_interfaces(self):
        """
        Test interface list XML file with no interfaces present.
        """
        path_output = self.base_output_dir.joinpath("interfaces_extended.json")

        # Run
        converter = ConvertInterfaceListToJSON(
            path_txt=str(
                self.base_input_dir.joinpath("interfaces_extended_no_interfaces.txt")
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        self.assertFalse(
            path_output.exists(),
            msg="Interface list JSON created when no interfaces were present.",
        )


class TestConvertAssemblyListToJSON(TestCase):
    """
    Tests for ConvertAssemblyListToJSON class.
    """

    def setUp(self):
        super().setUp()

        self.base_data_dir = Path("tests/data/")
        self.base_expected_dir = self.base_data_dir.joinpath(
            "expected_output/list_results/"
        )
        self.base_output_dir = self.base_data_dir.joinpath(
            "actual_output/list_results/"
        )
        self.path_base_input = self.base_data_dir.joinpath("mock_data/list_results/")

        # Remove any existing output data
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        remove_files(self.base_output_dir)

    def test_parse_arbitrary_minimal(self):
        """
        Arbitrary minimal interface list XML file.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_minimal.json")
        converter = ConvertAssemblyListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "arbitrary_data_minimal", "assemblies_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_minimal", "assemblies_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary minimal interface list XML->JSON not parsed correctly.",
        )

    def test_parse_6nxr(self):
        """
        Test assembly list XML file from PISA analysis of 6nxr.
        """
        path_output = self.base_output_dir.joinpath("assemblies_extended.json")
        # Run
        converter = ConvertAssemblyListToJSON(
            path_txt=str(
                self.path_base_input.joinpath("6nxr", "assemblies_extended.txt")
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr", "assemblies_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="6nxr assembly list XML->JSON not parsed correctly.",
        )

    def test_parse_variable_column_widths(self):
        """
        Test assembly list XML file with variable column widths.
        """
        path_output = self.base_output_dir.joinpath(
            "assemblies_extended_variable_widths.json"
        )

        # Run
        converter = ConvertAssemblyListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "6nxr_variable_column_widths", "assemblies_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr_variable_column_width", "assemblies_extended.json"
        )

        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Variable column widths assembly list XML->JSON not parsed correctly.",
        )

    def test_parse_arbitrary_maximal(self):
        """
        Arbitrary maximal complexity assembly list XML file. Should cover all edge
        cases.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_maximal.json")
        converter = ConvertAssemblyListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "arbitrary_data_maximal", "assemblies_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_maximal", "assemblies_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary maximal assembly list XML->JSON not parsed correctly.",
        )


class TestConvertMonomerListToJSON(TestCase):
    """
    Tests for ConvertMonomerListToJSON class.
    """

    def setUp(self):
        super().setUp()

        self.maxDiff = None

        self.base_data_dir = Path("tests/data/")
        self.base_expected_dir = self.base_data_dir.joinpath(
            "expected_output/list_results/"
        )
        self.base_output_dir = self.base_data_dir.joinpath(
            "actual_output/list_results/"
        )
        self.path_base_input = self.base_data_dir.joinpath("mock_data/list_results/")

        # Remove any existing output data
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        remove_files(self.base_output_dir)

    def test_parse_arbitrary_minimal(self):
        """
        Arbitrary minimal monomer list XML file.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_minimal.json")
        converter = ConvertComponentsListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "arbitrary_data_minimal", "monomers_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_minimal", "monomers_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary minimal monomer list XML->JSON not parsed correctly.",
        )

    def test_parse_6nxr(self):
        """
        Test monomer list XML file from PISA analysis of 6nxr.
        """
        path_output = self.base_output_dir.joinpath("monomers_extended.json")
        # Run
        converter = ConvertComponentsListToJSON(
            path_txt=str(
                self.path_base_input.joinpath("6nxr", "monomers_extended.txt")
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr", "monomers_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="6nxr monomer list XML->JSON not parsed correctly.",
        )

    def test_parse_variable_column_widths(self):
        """
        Test monomer list XML file with variable column widths.
        """
        path_output = self.base_output_dir.joinpath(
            "monomers_extended_variable_widths.json"
        )

        # Run
        converter = ConvertComponentsListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "6nxr_variable_column_widths", "monomers_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "6nxr", "monomers_extended.json"
        )

        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Variable column widths monomer list XML->JSON not parsed correctly.",
        )

    def test_parse_arbitrary_maximal(self):
        """
        Arbitrary maximal complexity monomer list XML file. Should cover all edge
        cases.
        """
        # Run
        path_output = self.base_output_dir.joinpath("arbitrary_data_maximal.json")
        converter = ConvertComponentsListToJSON(
            path_txt=str(
                self.path_base_input.joinpath(
                    "arbitrary_data_maximal", "monomers_extended.txt"
                )
            ),
            path_json=str(path_output),
        )
        converter.parse()

        # Check
        path_expected = self.base_expected_dir.joinpath(
            "arbitrary_data_maximal", "monomers_extended.json"
        )
        json_expected = json.loads(path_expected.read_text().strip())
        json_actual = json.loads(path_output.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Arbitrary maximal monomer list XML->JSON not parsed correctly.",
        )


class TestCompileInterfaceSummaryJSON(TestCase):
    """
    Tests for CompileInterfaceSummaryJSON class.
    """

    def setUp(self):
        super().setUp()

        self.maxDiff = None

        self.base_input_dir = Path("tests/data/")

        # Remove any existing output data
        output_path = Path("tests/data/actual_output/")
        output_path.mkdir(parents=True, exist_ok=True)
        remove_files(output_path)

    def test_compile_interface_summary_json(self):
        """
        Test for multiple interfaces present.
        """

        self.input_interface_jsons = self.base_input_dir.joinpath(
            "mock_data",
            "interface_summary_parser",
            "interfaces_minified",
        )
        self.input_assembly_json = self.base_input_dir.joinpath(
            "mock_data",
            "interface_summary_parser",
            "assemblies.json",
        )
        self.output_json = self.base_input_dir.joinpath(
            "actual_output",
            "interface_summary.json",
        )
        self.expected_json = self.base_input_dir.joinpath(
            "expected_output",
            "interface_summary.json",
        )

        # Run
        self.compiler = CompileInterfaceSummaryJSON(
            path_interface_jsons=self.input_interface_jsons,
            path_assembly_json=self.input_assembly_json,
            path_output_json=self.output_json,
        )
        self.compiler.parse()

        # Check
        json_expected = json.loads(self.expected_json.read_text().strip())
        json_actual = json.loads(self.output_json.read_text().strip())

        self.assertDictEqual(
            json_expected,
            json_actual,
            msg="Interface summary JSON not compiled correctly.",
        )

    def test_compile_interface_summary_json_no_interfaces(self):
        """
        Test for no interfaces present. Should not create output JSON.
        """

        self.input_interface_jsons = self.base_input_dir.joinpath(
            "mock_data",
            "interface_summary_parser",
            "interfaces_none",
        )
        self.input_assembly_json = self.base_input_dir.joinpath(
            "mock_data",
            "interface_summary_parser",
            "assemblies.json",
        )
        self.output_json = self.base_input_dir.joinpath(
            "actual_output",
            "interface_summary_no_interfaces.json",
        )

        # Run
        self.compiler = CompileInterfaceSummaryJSON(
            path_interface_jsons=self.input_interface_jsons,
            path_assembly_json=self.input_assembly_json,
            path_output_json=self.output_json,
        )
        self.compiler.parse()

        # Check
        self.assertFalse(
            self.output_json.exists(),
            msg="Interface summary JSON created when no interfaces were present.",
        )
