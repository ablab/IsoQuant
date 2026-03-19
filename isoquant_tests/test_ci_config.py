############################################################################
# Copyright (c) 2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Unit tests for CI test config loading (YAML and TSV formats)."""

import os
import sys
import tempfile
import shutil
import unittest

import yaml

# Add isoquant_tests/github to path so we can import the modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "github"))

from run_pipeline import (
    load_tsv_config as pipeline_load_tsv,
    load_yaml_config as pipeline_load_yaml,
    load_config as pipeline_load_config,
    check_value,
)
from run_barcode_test import (
    load_config as barcode_load_config,
)
from cfg2yaml import (
    detect_barcode_mode,
    convert_cfg_to_yaml,
    load_baseline_file,
)
from update_defaults import (
    update_yaml_baselines,
    load_tsv_config as defaults_load_tsv,
)


class TestLoadTsvConfig(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_basic_loading(self):
        path = os.path.join(self.test_dir, "test.cfg")
        with open(path, "w") as f:
            f.write("name\tmy_test\n")
            f.write("output\t/tmp/out\n")
            f.write("datatype\tont\n")
        result = pipeline_load_tsv(path)
        self.assertEqual(result["name"], "my_test")
        self.assertEqual(result["output"], "/tmp/out")
        self.assertEqual(result["datatype"], "ont")

    def test_comments_skipped(self):
        path = os.path.join(self.test_dir, "test.cfg")
        with open(path, "w") as f:
            f.write("# this is a comment\n")
            f.write("name\tmy_test\n")
            f.write("# another comment\n")
        result = pipeline_load_tsv(path)
        self.assertEqual(len(result), 1)
        self.assertEqual(result["name"], "my_test")

    def test_blank_lines_skipped(self):
        path = os.path.join(self.test_dir, "test.cfg")
        with open(path, "w") as f:
            f.write("name\tmy_test\n")
            f.write("\n")
            f.write("output\t/tmp\n")
        result = pipeline_load_tsv(path)
        self.assertEqual(len(result), 2)

    def test_missing_file_returns_empty(self):
        result = pipeline_load_tsv("/nonexistent/path.cfg")
        self.assertEqual(result, {})


class TestLoadYamlConfig(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write_yaml(self, data: dict, name: str = "test.yaml") -> str:
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        return path

    def test_basic_loading(self):
        path = self._write_yaml({
            "name": "my_test",
            "output": "/tmp/out",
            "datatype": "ont",
        })
        config, baselines = pipeline_load_yaml(path)
        self.assertEqual(config["name"], "my_test")
        self.assertEqual(config["output"], "/tmp/out")
        self.assertEqual(baselines, {})

    def test_baselines_extracted(self):
        path = self._write_yaml({
            "name": "my_test",
            "baselines": {
                "transcripts": {"full_recall": 84.7, "full_precision": 96.6},
                "performance": {"max_rss": 1234.0},
            },
        })
        config, baselines = pipeline_load_yaml(path)
        self.assertNotIn("baselines", config)
        self.assertIn("transcripts", baselines)
        self.assertAlmostEqual(baselines["transcripts"]["full_recall"], 84.7)
        self.assertIn("performance", baselines)

    def test_values_converted_to_strings(self):
        path = self._write_yaml({
            "name": "test",
            "threads": 8,
            "tolerance": 0.01,
        })
        config, _ = pipeline_load_yaml(path)
        self.assertEqual(config["threads"], "8")
        self.assertEqual(config["tolerance"], "0.01")

    def test_no_baselines_key(self):
        path = self._write_yaml({"name": "test"})
        config, baselines = pipeline_load_yaml(path)
        self.assertEqual(baselines, {})

    def test_missing_file_returns_empty(self):
        config, baselines = pipeline_load_yaml("/nonexistent/path.yaml")
        self.assertEqual(config, {})
        self.assertEqual(baselines, {})

    def test_empty_file(self):
        path = os.path.join(self.test_dir, "empty.yaml")
        with open(path, "w") as f:
            f.write("")
        config, baselines = pipeline_load_yaml(path)
        self.assertEqual(config, {})
        self.assertEqual(baselines, {})


class TestLoadConfigDispatch(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_yaml_extension(self):
        path = os.path.join(self.test_dir, "test.yaml")
        with open(path, "w") as f:
            yaml.dump({"name": "yaml_test", "baselines": {"transcripts": {"r": 1.0}}}, f)
        config, baselines = pipeline_load_config(path)
        self.assertEqual(config["name"], "yaml_test")
        self.assertIn("transcripts", baselines)

    def test_yml_extension(self):
        path = os.path.join(self.test_dir, "test.yml")
        with open(path, "w") as f:
            yaml.dump({"name": "yml_test"}, f)
        config, baselines = pipeline_load_config(path)
        self.assertEqual(config["name"], "yml_test")
        self.assertEqual(baselines, {})

    def test_cfg_extension(self):
        path = os.path.join(self.test_dir, "test.cfg")
        with open(path, "w") as f:
            f.write("name\tcfg_test\n")
        config, baselines = pipeline_load_config(path)
        self.assertEqual(config["name"], "cfg_test")
        self.assertEqual(baselines, {})

    def test_barcode_load_config_yaml(self):
        path = os.path.join(self.test_dir, "bc.yaml")
        with open(path, "w") as f:
            yaml.dump({
                "name": "bc_test",
                "mode": "custom_sc",
                "baselines": {"barcode": {"precision": 99.5}},
            }, f)
        config, baselines = barcode_load_config(path)
        self.assertEqual(config["name"], "bc_test")
        self.assertAlmostEqual(baselines["barcode"]["precision"], 99.5)


class TestCheckValue(unittest.TestCase):
    def test_within_tolerance(self):
        self.assertEqual(check_value(100.0, 100.5, "metric"), 0)

    def test_below_tolerance(self):
        self.assertNotEqual(check_value(100.0, 90.0, "metric"), 0)

    def test_above_tolerance(self):
        self.assertNotEqual(check_value(100.0, 110.0, "metric"), 0)

    def test_exact_match(self):
        self.assertEqual(check_value(50.0, 50.0, "metric"), 0)

    def test_negative_exact_match(self):
        self.assertEqual(check_value(-1.0, -1.0, "metric"), 0)

    def test_negative_mismatch(self):
        self.assertNotEqual(check_value(-1.0, -2.0, "metric"), 0)

    def test_custom_tolerance(self):
        # 20% tolerance: 100 * 0.8 = 80, so 85 should pass
        self.assertEqual(check_value(100.0, 85.0, "metric", 0.2), 0)
        # but 70 should fail
        self.assertNotEqual(check_value(100.0, 70.0, "metric", 0.2), 0)


class TestDetectBarcodeMode(unittest.TestCase):
    def test_pipeline_config(self):
        config = {"genome": "/path/fa", "genedb": "/path/gtf", "datatype": "ont"}
        self.assertFalse(detect_barcode_mode(config))

    def test_barcode_config(self):
        config = {"mode": "custom_sc", "molecule": "10x.mdf", "qa_mode": "tenX_v3"}
        self.assertTrue(detect_barcode_mode(config))

    def test_sc_pipeline_config(self):
        # SC pipeline has both barcode-like isoquant_options and pipeline keys
        # but barcode keys are inside isoquant_options, not top-level
        config = {"genome": "/path/fa", "genedb": "/path/gtf", "datatype": "ont",
                  "isoquant_options": '"-m tenX_v3"'}
        self.assertFalse(detect_barcode_mode(config))

    def test_mixed_keys_pipeline_wins(self):
        # If both barcode and pipeline keys are present, it's a pipeline config
        config = {"genome": "/path/fa", "genedb": "/path/gtf", "mode": "custom_sc"}
        self.assertFalse(detect_barcode_mode(config))


class TestCfg2Yaml(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_pipeline_conversion(self):
        # Create cfg
        cfg_path = os.path.join(self.test_dir, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\tmy_pipeline\n")
            f.write("genome\t/path/to/genome.fa\n")
            f.write("genedb\t/path/to/genes.gtf\n")
            f.write("datatype\tont\n")
            f.write("etalon\ttest.etl\n")

        # Create etalon file
        etl_path = os.path.join(self.test_dir, "test.etl")
        with open(etl_path, "w") as f:
            f.write("full_recall\t84.70\n")
            f.write("full_precision\t96.60\n")

        yaml_path, yaml_data = convert_cfg_to_yaml(cfg_path)
        self.assertEqual(yaml_path, os.path.join(self.test_dir, "test.yaml"))
        self.assertEqual(yaml_data["name"], "my_pipeline")
        self.assertNotIn("etalon", yaml_data)
        self.assertIn("baselines", yaml_data)
        self.assertAlmostEqual(yaml_data["baselines"]["transcripts"]["full_recall"], 84.7)
        self.assertAlmostEqual(yaml_data["baselines"]["transcripts"]["full_precision"], 96.6)

    def test_barcode_conversion(self):
        cfg_path = os.path.join(self.test_dir, "bc.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\tmy_barcode\n")
            f.write("mode\tcustom_sc\n")
            f.write("molecule\t10x.mdf\n")
            f.write("etalon\tbc.etl\n")

        etl_path = os.path.join(self.test_dir, "bc.etl")
        with open(etl_path, "w") as f:
            f.write("precision\t99.89\n")
            f.write("recall\t86.82\n")

        yaml_path, yaml_data = convert_cfg_to_yaml(cfg_path)
        self.assertNotIn("etalon", yaml_data)
        self.assertIn("baselines", yaml_data)
        self.assertIn("barcode", yaml_data["baselines"])
        self.assertAlmostEqual(yaml_data["baselines"]["barcode"]["precision"], 99.89)

    def test_force_barcode_mode(self):
        cfg_path = os.path.join(self.test_dir, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\ttest\n")
            f.write("etalon\ttest.etl\n")

        etl_path = os.path.join(self.test_dir, "test.etl")
        with open(etl_path, "w") as f:
            f.write("metric1\t1.0\n")

        _, yaml_data = convert_cfg_to_yaml(cfg_path, force_barcode=True)
        self.assertIn("barcode", yaml_data["baselines"])
        self.assertNotIn("transcripts", yaml_data.get("baselines", {}))

    def test_missing_etalon_no_baselines(self):
        cfg_path = os.path.join(self.test_dir, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\ttest\n")
            f.write("genome\t/path/genome.fa\n")
            f.write("datatype\tont\n")

        _, yaml_data = convert_cfg_to_yaml(cfg_path)
        self.assertNotIn("baselines", yaml_data)


class TestLoadBaselineFile(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_basic(self):
        path = os.path.join(self.test_dir, "test.etl")
        with open(path, "w") as f:
            f.write("recall\t84.70\n")
            f.write("precision\t96.60\n")
        result = load_baseline_file(path)
        self.assertAlmostEqual(result["recall"], 84.7)
        self.assertAlmostEqual(result["precision"], 96.6)

    def test_missing_file(self):
        result = load_baseline_file("/nonexistent/file.etl")
        self.assertEqual(result, {})

    def test_comments_and_blanks(self):
        path = os.path.join(self.test_dir, "test.etl")
        with open(path, "w") as f:
            f.write("# header\n")
            f.write("metric1\t1.0\n")
            f.write("\n")
            f.write("metric2\t2.0\n")
        result = load_baseline_file(path)
        self.assertEqual(len(result), 2)


class TestUpdateYamlBaselines(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.run_name = "my_test"
        # Create pipeline output structure
        self.pipeline_dir = os.path.join(self.test_dir, "output", self.run_name)
        os.makedirs(os.path.join(self.pipeline_dir, "gffcompare"))

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def _write_yaml_config(self, data: dict) -> str:
        path = os.path.join(self.test_dir, "config.yaml")
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        return path

    def _write_tsv(self, path: str, data: dict) -> None:
        with open(path, "w") as f:
            for k, v in data.items():
                f.write("%s\t%.2f\n" % (k, v))

    def test_pipeline_baselines_updated(self):
        config_path = self._write_yaml_config({
            "name": self.run_name,
            "baselines": {
                "transcripts": {"full_recall": 80.0, "full_precision": 90.0},
            },
        })
        # Write new values
        self._write_tsv(
            os.path.join(self.pipeline_dir, "gffcompare", "new_gtf_etalon.tsv"),
            {"full_recall": 85.0, "full_precision": 95.0},
        )

        base_output_dir = os.path.join(self.test_dir, "output")
        update_yaml_baselines(config_path, self.run_name, base_output_dir)

        with open(config_path) as f:
            updated = yaml.safe_load(f)
        self.assertAlmostEqual(updated["baselines"]["transcripts"]["full_recall"], 85.0)
        self.assertAlmostEqual(updated["baselines"]["transcripts"]["full_precision"], 95.0)

    def test_barcode_baselines_updated(self):
        config_path = self._write_yaml_config({
            "name": self.run_name,
            "baselines": {
                "barcode": {"precision": 99.0, "recall": 85.0},
            },
        })
        # Barcode etalon lives at {base_output_dir}/{run_name}.new_etalon.tsv
        base_output_dir = os.path.join(self.test_dir, "output")
        self._write_tsv(
            os.path.join(base_output_dir, self.run_name + ".new_etalon.tsv"),
            {"precision": 99.5, "recall": 87.0},
        )

        update_yaml_baselines(config_path, self.run_name, base_output_dir)

        with open(config_path) as f:
            updated = yaml.safe_load(f)
        self.assertAlmostEqual(updated["baselines"]["barcode"]["precision"], 99.5)
        self.assertAlmostEqual(updated["baselines"]["barcode"]["recall"], 87.0)

    def test_missing_output_skipped(self):
        config_path = self._write_yaml_config({
            "name": self.run_name,
            "baselines": {
                "transcripts": {"full_recall": 80.0},
                "performance": {"max_rss": 1000.0},
            },
        })
        # Only write transcripts output, not performance
        self._write_tsv(
            os.path.join(self.pipeline_dir, "gffcompare", "new_gtf_etalon.tsv"),
            {"full_recall": 85.0},
        )

        base_output_dir = os.path.join(self.test_dir, "output")
        update_yaml_baselines(config_path, self.run_name, base_output_dir)

        with open(config_path) as f:
            updated = yaml.safe_load(f)
        # transcripts updated
        self.assertAlmostEqual(updated["baselines"]["transcripts"]["full_recall"], 85.0)
        # performance unchanged
        self.assertAlmostEqual(updated["baselines"]["performance"]["max_rss"], 1000.0)

    def test_non_baseline_config_preserved(self):
        config_path = self._write_yaml_config({
            "name": self.run_name,
            "output": "/some/path",
            "datatype": "ont",
            "baselines": {
                "transcripts": {"full_recall": 80.0},
            },
        })
        self._write_tsv(
            os.path.join(self.pipeline_dir, "gffcompare", "new_gtf_etalon.tsv"),
            {"full_recall": 85.0},
        )

        base_output_dir = os.path.join(self.test_dir, "output")
        update_yaml_baselines(config_path, self.run_name, base_output_dir)

        with open(config_path) as f:
            updated = yaml.safe_load(f)
        self.assertEqual(updated["name"], self.run_name)
        self.assertEqual(updated["output"], "/some/path")
        self.assertEqual(updated["datatype"], "ont")


class TestYamlTsvRoundTrip(unittest.TestCase):
    """Verify that YAML and TSV configs produce equivalent results when loaded."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_config_values_match(self):
        # Create a TSV config
        cfg_path = os.path.join(self.test_dir, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\tmy_test\n")
            f.write("output\t/tmp/output\n")
            f.write("run_type\ttranscripts,quantification_known\n")
            f.write("genome\t/path/to/genome.fa\n")
            f.write("datatype\tont\n")
            f.write('isoquant_options\t"-t 8 --complete_genedb --force"\n')

        # Create equivalent YAML config
        yaml_path = os.path.join(self.test_dir, "test.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump({
                "name": "my_test",
                "output": "/tmp/output",
                "run_type": "transcripts,quantification_known",
                "genome": "/path/to/genome.fa",
                "datatype": "ont",
                "isoquant_options": '"-t 8 --complete_genedb --force"',
            }, f, default_flow_style=False, sort_keys=False)

        cfg_config, cfg_baselines = pipeline_load_config(cfg_path)
        yaml_config, yaml_baselines = pipeline_load_config(yaml_path)

        for key in cfg_config:
            self.assertEqual(cfg_config[key], yaml_config.get(key),
                             "Mismatch for key %s" % key)

    def test_isoquant_options_processing(self):
        """Verify isoquant_options string processes identically after loading."""
        cfg_path = os.path.join(self.test_dir, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write("name\ttest\n")
            f.write('isoquant_options\t"-t 8 --complete_genedb --force"\n')

        yaml_path = os.path.join(self.test_dir, "test.yaml")
        with open(yaml_path, "w") as f:
            yaml.dump({
                "name": "test",
                "isoquant_options": '"-t 8 --complete_genedb --force"',
            }, f)

        cfg_config, _ = pipeline_load_config(cfg_path)
        yaml_config, _ = pipeline_load_config(yaml_path)

        cfg_opts = cfg_config["isoquant_options"].replace('"', '').split()
        yaml_opts = yaml_config["isoquant_options"].replace('"', '').split()
        self.assertEqual(cfg_opts, yaml_opts)


if __name__ == "__main__":
    unittest.main()
