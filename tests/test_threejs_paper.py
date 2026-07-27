import re
import shutil
import subprocess
import unittest
from pathlib import Path

from oviz.paper import Paper
from oviz.threejs_figure import ThreeJSFigure
from oviz.threejs_paper import normalize_paper_spec, resolve_paper_state_bindings


class PaperSpecTests(unittest.TestCase):
    def test_normalize_defaults_and_clamps(self):
        spec = normalize_paper_spec(None)
        self.assertEqual(spec["schema_version"], 1)
        self.assertFalse(spec["available"])
        self.assertFalse(spec["enabled"])
        self.assertEqual(spec["sections"], [])

        spec = normalize_paper_spec(
            {
                "title": "T",
                "panel": {"width_fraction": 0.95},
                "sync": {"mode": "bogus", "reading_line": 5.0, "resume_policy": "nah"},
                "sections": [
                    {
                        "title": "Intro",
                        "blocks": [
                            {"type": "html", "html": "<p>x</p>", "anchor": {"state": "Overview"}},
                            {"type": "figure", "image_data_url": "data:image/png;base64,QUJD"},
                            {"type": "weird"},
                        ],
                    }
                ],
            }
        )
        self.assertTrue(spec["available"])
        self.assertTrue(spec["enabled"])
        self.assertEqual(spec["panel"]["width_fraction"], 0.60)
        self.assertEqual(spec["sync"]["mode"], "scroll")
        self.assertEqual(spec["sync"]["resume_policy"], "next-anchor")
        blocks = spec["sections"][0]["blocks"]
        self.assertEqual(blocks[0]["anchor"]["state"], "Overview")
        self.assertIsNone(blocks[0]["anchor"]["state_id"])
        self.assertEqual(blocks[1]["type"], "figure")
        self.assertEqual(blocks[2]["type"], "html")

    def test_large_figures_are_content_addressed(self):
        big = "data:image/png;base64," + "A" * 5000
        spec = normalize_paper_spec(
            {
                "sections": [
                    {"blocks": [
                        {"type": "figure", "image_data_url": big},
                        {"type": "figure", "image_data_url": big},
                    ]}
                ]
            }
        )
        blocks = spec["sections"][0]["blocks"]
        self.assertIsInstance(blocks[0]["image_data_url"], dict)
        self.assertIn("__oviz_asset_ref__", blocks[0]["image_data_url"])
        self.assertEqual(blocks[0]["image_data_url"], blocks[1]["image_data_url"])
        self.assertEqual(len(spec["assets"]), 1)
        # Normalization must be idempotent: ThreeJSFigure re-normalizes the
        # paper spec, and a second pass must keep asset references as objects
        # rather than coercing them to strings.
        again = normalize_paper_spec(spec)
        again_blocks = again["sections"][0]["blocks"]
        self.assertEqual(again_blocks[0]["image_data_url"], blocks[0]["image_data_url"])
        self.assertEqual(again["assets"], spec["assets"])

    def test_state_binding_resolution_by_name(self):
        states = {"items": [{"id": "state-a", "name": "Overview"}]}
        spec = normalize_paper_spec(
            {"sections": [{"blocks": [{"html": "<p>x</p>", "anchor": {"state": "overview"}}]}]}
        )
        resolved = resolve_paper_state_bindings(spec, states)
        anchor = resolved["sections"][0]["blocks"][0]["anchor"]
        self.assertEqual(anchor["state_id"], "state-a")

        with self.assertRaises(ValueError) as raised:
            resolve_paper_state_bindings(
                normalize_paper_spec(
                    {"sections": [{"blocks": [{"html": "x", "anchor": {"state": "Missing"}}]}]}
                ),
                states,
            )
        self.assertIn("Missing", str(raised.exception))
        self.assertIn("Overview", str(raised.exception))

    def test_paper_builder_round_trip(self):
        paper = Paper("Title", authors=["A. Author"], link_url="https://example.org")
        paper.add_section("Intro")
        paper.add_paragraph("Hello.", state="Overview", label="Start")
        spec = paper.to_spec(states={"items": [{"id": "s1", "name": "Overview"}]})
        self.assertTrue(spec["available"])
        anchor = spec["sections"][0]["blocks"][0]["anchor"]
        self.assertEqual(anchor["state_id"], "s1")


class PaperRuntimeTests(unittest.TestCase):
    def _html(self, paper=None):
        scene = {"width": 640, "height": 480, "frames": [], "initial_state": {}}
        if paper is not None:
            scene["paper"] = paper
        return ThreeJSFigure(scene).to_html(compress_scene_spec=False)

    def test_runtime_and_panel_ship_when_paper_available(self):
        paper = Paper("T", enabled=True)
        paper.add_section("S")
        paper.add_paragraph("text")
        html = self._html(paper.to_spec())
        self.assertIn("function initializeOvizPaper()", html)
        self.assertIn('class="oviz-three-paper-panel"', html)
        self.assertIn("Paper ▸", html)
        self.assertIn("function ovizPaperQueueStateJump(target)", html)
        self.assertIn("ovizPaperPauseFromSceneInteraction", html)
        # the paper panel participates in the interaction guards
        self.assertIn('closest(".oviz-three-paper-panel")', html)

    def test_paper_button_absent_without_paper(self):
        html = self._html(None)
        self.assertNotIn("Paper ▸", html)
        self.assertIn('class="oviz-three-paper-panel"', html)  # markup ships, runtime removes
        # The collapse/expand edge tab must never flash or linger in figures
        # without a paper: hidden in markup, removed by the runtime.
        self.assertIn('class="oviz-three-paper-expand-tab" type="button" hidden', html)
        self.assertIn('strayExpandTab.remove()', html)


class PaperDemoArtifactTests(unittest.TestCase):
    ARTIFACT = Path(__file__).with_name("main_figure_paper_demo.html")

    def test_demo_artifact_content(self):
        if not self.ARTIFACT.exists():
            self.skipTest("paper demo artifact not built")
        html = self.ARTIFACT.read_text(encoding="utf-8")
        self.assertLess(self.ARTIFACT.stat().st_size, 100 * 1024 * 1024)
        self.assertIn("Paper ▸", html)
        self.assertIn("function initializeOvizPaper()", html)
        if shutil.which("node"):
            scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
            viewer_script = next(
                script for script in scripts if "/*__SCENE_SPEC_START__*/" in script
            )
            script_path = self.ARTIFACT.with_name("paper-demo-check.js")
            try:
                script_path.write_text(viewer_script, encoding="utf-8")
                subprocess.run(["node", "--check", str(script_path)], check=True)
            finally:
                script_path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
